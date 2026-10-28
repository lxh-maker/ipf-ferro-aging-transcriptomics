# -*- coding: utf-8 -*-
options(encoding = "UTF-8", stringsAsFactors = FALSE)
try(Sys.setlocale(category = "LC_CTYPE", locale = "Chinese_China.65001"), silent = TRUE)

suppressPackageStartupMessages(library(e1071))

workDir <- Sys.getenv("IPF_BULK_WORKDIR", unset = getwd())
workDir <- normalizePath(workDir, winslash = "/", mustWork = TRUE)
setwd(workDir)

inputFile <- file.path(workDir, "Sample Type Matrix.csv")
geneListFile <- file.path(workDir, "Intersection_DEG_ferro-aging.csv")
outer_folds <- 5L
outer_repeats <- 20L
base_seed <- 12345L

strip_bom <- function(x) sub("^\ufeff", "", x)
clean_gene_vector <- function(x) {
  x <- strip_bom(trimws(as.character(x)))
  x <- gsub("^['\"]|['\"]$", "", x)
  unique(x[!is.na(x) & x != ""])
}

read_gene_list <- function(file_path) {
  if (!file.exists(file_path)) stop("基因列表文件不存在：", file_path)
  gene_df <- read.csv(file_path, header = TRUE, check.names = FALSE,
                      stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  names(gene_df) <- strip_bom(names(gene_df))
  gene_col <- which(tolower(gsub("[^A-Za-z0-9]+", "", names(gene_df))) %in%
                      c("gene", "genes", "genename", "genesymbol", "symbol"))
  if (length(gene_col) == 0) gene_col <- 1L
  genes <- clean_gene_vector(gene_df[[gene_col[1]]])
  if (length(genes) == 0) stop("基因列表文件中没有有效基因。")
  genes
}

order_binary_group <- function(group) {
  group <- factor(group)
  if (nlevels(group) != 2) stop("SVM-RFE仅支持二分类；当前分组：", paste(levels(group), collapse = ", "))
  lev <- levels(group)
  control_hit <- grepl("^(con|control|normal|healthy)$", lev, ignore.case = TRUE)
  if (sum(control_hit) == 1) lev <- c(lev[control_hit], lev[!control_hit])
  factor(as.character(group), levels = lev)
}

make_stratified_folds <- function(y, k = 5, seed = 1L) {
  y <- factor(y)
  k <- min(k, min(table(y)))
  if (k < 3) stop("每组样本数太少，无法进行至少3折交叉验证。")
  set.seed(seed)
  foldid <- integer(length(y))
  for (lev in levels(y)) {
    idx <- which(y == lev)
    foldid[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  foldid
}

fit_scaler <- function(x) {
  center <- colMeans(x, na.rm = TRUE)
  scale_value <- apply(x, 2, sd, na.rm = TRUE)
  scale_value[!is.finite(scale_value) | scale_value == 0] <- 1
  list(center = center, scale = scale_value)
}

apply_scaler <- function(x, scaler) {
  out <- sweep(x, 2, scaler$center, "-")
  sweep(out, 2, scaler$scale, "/")
}

class_weights <- function(y) {
  tab <- table(y)
  setNames(length(y) / (2 * as.numeric(tab)), names(tab))
}

# 仅在训练折中完成递归特征消除，避免测试折信息泄漏。
linear_svm_rfe_rank <- function(x_train, y_train) {
  surviving <- seq_len(ncol(x_train))
  eliminated <- integer(0)
  cw <- class_weights(y_train)
  while (length(surviving) > 0) {
    if (length(surviving) == 1) {
      eliminated <- c(eliminated, surviving)
      break
    }
    fit <- svm(
      x = x_train[, surviving, drop = FALSE], y = y_train,
      type = "C-classification", kernel = "linear", cost = 1,
      scale = FALSE, class.weights = cw
    )
    weight <- drop(t(fit$coefs) %*% fit$SV)
    worst_local <- which.min(weight^2)
    eliminated <- c(eliminated, surviving[worst_local])
    surviving <- surviving[-worst_local]
  }
  rev(eliminated)
}

calc_metrics <- function(truth, pred, positive) {
  truth <- factor(truth, levels = levels(truth))
  pred <- factor(pred, levels = levels(truth))
  negative <- setdiff(levels(truth), positive)
  sensitivity <- mean(pred[truth == positive] == positive)
  specificity <- mean(pred[truth == negative] == negative)
  accuracy <- mean(pred == truth)
  c(Accuracy = accuracy, Sensitivity = sensitivity, Specificity = specificity,
    BalancedAccuracy = mean(c(sensitivity, specificity)))
}

if (!file.exists(inputFile)) stop("表达数据文件不存在：", inputFile)
geneList <- read_gene_list(geneListFile)
raw <- read.csv(inputFile, header = TRUE, check.names = FALSE, row.names = 1,
                fileEncoding = "UTF-8-BOM")
availableGenes <- intersect(geneList, rownames(raw))
if (length(availableGenes) < 2) stop("表达矩阵中匹配到的候选基因少于2个。")

x <- t(as.matrix(raw[availableGenes, , drop = FALSE]))
storage.mode(x) <- "double"
if (anyNA(x)) stop("候选基因表达矩阵中存在NA，请先处理缺失值。")
y <- order_binary_group(sub("^.*_", "", rownames(x)))
positive <- levels(y)[2]
n <- nrow(x)
p <- ncol(x)
outer_folds <- min(outer_folds, min(table(y)))

message(sprintf("样本=%d，基因=%d，分组=%s", n, p,
                paste(names(table(y)), as.integer(table(y)), sep = "=", collapse = "; ")))
message(sprintf("开始%d次重复、%d折分层交叉验证...", outer_repeats, outer_folds))

metric_array <- array(NA_real_, dim = c(outer_repeats, p, 4),
                      dimnames = list(NULL, paste0("N", seq_len(p)),
                                      c("Accuracy", "Sensitivity", "Specificity", "BalancedAccuracy")))
rank_records <- vector("list", outer_repeats * outer_folds)
pred_array <- array(NA_integer_, dim = c(n, p, outer_repeats),
                    dimnames = list(rownames(x), paste0("N", seq_len(p)), paste0("Repeat", seq_len(outer_repeats))))
fold_records <- list()
record_id <- 0L

for (r in seq_len(outer_repeats)) {
  foldid <- make_stratified_folds(y, outer_folds, base_seed + r)
  fold_records[[r]] <- data.frame(Sample = rownames(x), Group = as.character(y),
                                  Repeat = r, Fold = foldid)
  for (fold in seq_len(outer_folds)) {
    record_id <- record_id + 1L
    test_idx <- which(foldid == fold)
    train_idx <- which(foldid != fold)
    scaler <- fit_scaler(x[train_idx, , drop = FALSE])
    x_train <- apply_scaler(x[train_idx, , drop = FALSE], scaler)
    x_test <- apply_scaler(x[test_idx, , drop = FALSE], scaler)
    y_train <- y[train_idx]
    ranking <- linear_svm_rfe_rank(x_train, y_train)
    rank_position <- integer(p)
    rank_position[ranking] <- seq_len(p)
    rank_records[[record_id]] <- rank_position
    cw <- class_weights(y_train)

    for (feature_n in seq_len(p)) {
      use_idx <- ranking[seq_len(feature_n)]
      fit <- svm(
        x = x_train[, use_idx, drop = FALSE], y = y_train,
        type = "C-classification", kernel = "linear", cost = 1,
        scale = FALSE, class.weights = cw
      )
      pred <- predict(fit, x_test[, use_idx, drop = FALSE])
      pred_array[test_idx, feature_n, r] <- as.integer(pred == positive)
    }
  }

  for (feature_n in seq_len(p)) {
    pred_factor <- factor(ifelse(pred_array[, feature_n, r] == 1, positive, levels(y)[1]),
                          levels = levels(y))
    metric_array[r, feature_n, ] <- calc_metrics(y, pred_factor, positive)
  }
  message(sprintf("已完成 %d/%d 次重复", r, outer_repeats))
}

rank_matrix <- do.call(rbind, rank_records)
colnames(rank_matrix) <- colnames(x)
avg_rank <- colMeans(rank_matrix)
sd_rank <- apply(rank_matrix, 2, sd)

metric_mean <- apply(metric_array, c(2, 3), mean, na.rm = TRUE)
metric_sd <- apply(metric_array, c(2, 3), sd, na.rm = TRUE)
metric_se <- metric_sd / sqrt(outer_repeats)
performance <- data.frame(
  Features = seq_len(p),
  Accuracy = metric_mean[, "Accuracy"],
  Sensitivity = metric_mean[, "Sensitivity"],
  Specificity = metric_mean[, "Specificity"],
  BalancedAccuracy = metric_mean[, "BalancedAccuracy"],
  BalancedAccuracy_SD = metric_sd[, "BalancedAccuracy"],
  BalancedAccuracy_SE = metric_se[, "BalancedAccuracy"]
)

# 1-SE规则：在接近最佳性能的候选中选择基因数最少的模型。
best_index <- which.max(performance$BalancedAccuracy)
one_se_cutoff <- performance$BalancedAccuracy[best_index] - performance$BalancedAccuracy_SE[best_index]
optimalFeatureCount <- min(performance$Features[performance$BalancedAccuracy >= one_se_cutoff])
optimalError <- 1 - performance$BalancedAccuracy[optimalFeatureCount]
optimalAccuracy <- performance$BalancedAccuracy[optimalFeatureCount]

final_order <- order(avg_rank, sd_rank)
selection_frequency <- colMeans(rank_matrix <= optimalFeatureCount)
topFeatures <- data.frame(
  FeatureName = colnames(x)[final_order],
  AvgRank = avg_rank[final_order],
  RankSD = sd_rank[final_order],
  SelectionFrequency = selection_frequency[final_order]
)
featureGenes <- topFeatures$FeatureName[seq_len(optimalFeatureCount)]

write.table(topFeatures, "feature_svm.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(featureGenes, "SVM-RFE.gene.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
write.csv(performance, "SVM_CV_Performance.csv", row.names = FALSE)
write.csv(do.call(rbind, fold_records), "SVM_Stratified_Folds.csv", row.names = FALSE)

oof_vote <- rowMeans(pred_array[, optimalFeatureCount, , drop = FALSE][, 1, ], na.rm = TRUE)
oof_class <- ifelse(oof_vote >= 0.5, positive, levels(y)[1])
write.csv(data.frame(Sample = rownames(x), Group = as.character(y),
                     OOF_DiseaseVote = oof_vote, OOF_Class = oof_class),
          "SVM_OOF_Predictions.csv", row.names = FALSE)

errors <- 1 - performance$BalancedAccuracy
pdf("errors.pdf", width = 5.5, height = 5)
plot(performance$Features, errors, type = "b", pch = 19, lwd = 2,
     xlab = "Number of Features", ylab = "Repeated CV Balanced Error",
     main = "SVM-RFE Generalization Error", ylim = range(c(errors, 0.5)))
abline(h = 0.5, lty = 2, col = "gray60")
points(optimalFeatureCount, optimalError, pch = 19, cex = 1.4, col = "red")
dev.off()

pdf("accuracy.pdf", width = 5.5, height = 5)
plot(performance$Features, performance$BalancedAccuracy, type = "b", pch = 19, lwd = 2,
     xlab = "Number of Features", ylab = "Repeated CV Balanced Accuracy",
     main = "SVM-RFE Cross-validated Performance", ylim = c(0, 1))
arrows(performance$Features,
       pmax(0, performance$BalancedAccuracy - performance$BalancedAccuracy_SE),
       performance$Features,
       pmin(1, performance$BalancedAccuracy + performance$BalancedAccuracy_SE),
       angle = 90, code = 3, length = 0.04, col = "gray50")
points(optimalFeatureCount, optimalAccuracy, pch = 19, cex = 1.4, col = "red")
dev.off()

pdf("combined_plots.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
plot(performance$Features, errors, type = "b", pch = 19, lwd = 2,
     xlab = "Number of Features", ylab = "Balanced Error", main = "Repeated Stratified CV")
points(optimalFeatureCount, optimalError, pch = 19, cex = 1.4, col = "red")
plot(performance$Features, performance$BalancedAccuracy, type = "b", pch = 19, lwd = 2,
     xlab = "Number of Features", ylab = "Balanced Accuracy", main = "Repeated Stratified CV",
     ylim = c(0, 1))
points(optimalFeatureCount, optimalAccuracy, pch = 19, cex = 1.4, col = "red")
dev.off()

summary_df <- data.frame(
  Item = c("samples", "groups", "candidate_genes", "outer_folds", "outer_repeats",
           "metric", "selection_rule", "selected_gene_count", "balanced_accuracy", "balanced_error"),
  Value = c(n, paste(names(table(y)), as.integer(table(y)), sep = "=", collapse = "; "), p,
            outer_folds, outer_repeats, "Balanced accuracy", "1-SE rule",
            optimalFeatureCount, optimalAccuracy, optimalError)
)
write.csv(summary_df, "SVM_Run_Summary.csv", row.names = FALSE)

cat("========================================\n")
cat("SVM-RFE分析完成（重复分层交叉验证）\n")
cat("最优特征数：", optimalFeatureCount, "\n", sep = "")
cat("平衡准确率：", sprintf("%.4f", optimalAccuracy), "\n", sep = "")
cat("筛选基因：", paste(featureGenes, collapse = ", "), "\n", sep = "")
cat("========================================\n")
