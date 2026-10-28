# -*- coding: utf-8 -*-
options(encoding = "UTF-8", stringsAsFactors = FALSE)
try(Sys.setlocale(category = "LC_CTYPE", locale = "Chinese_China.65001"), silent = TRUE)

suppressPackageStartupMessages({
  library(glmnet)
  library(pROC)
  library(ggplot2)
})

working_directory <- Sys.getenv("IPF_BULK_WORKDIR", unset = getwd())
working_directory <- normalizePath(working_directory, winslash = "/", mustWork = TRUE)
setwd(working_directory)

result_dir <- file.path(working_directory, "Analysis_Result")
if (!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)
input_expr <- file.path(working_directory, "Sample Type Matrix.csv")
gene_file <- file.path(working_directory, "Intersect_SVM_LASSO_RF_Genes.txt")
outer_folds <- 5L
outer_repeats <- 20L
base_seed <- 20261026L

strip_bom <- function(x) sub("^\ufeff", "", x)
clean_gene_vector <- function(x) {
  x <- strip_bom(trimws(as.character(x)))
  x <- gsub("^['\"]|['\"]$", "", x)
  x <- x[!is.na(x) & x != ""]
  header <- tolower(gsub("[^A-Za-z0-9]+", "", x)) %in% c("gene", "genes", "symbol", "genesymbol")
  unique(x[!header])
}

read_gene_list <- function(path) {
  if (!file.exists(path)) stop("未找到诊断基因文件：", path)
  genes <- clean_gene_vector(readLines(path, encoding = "UTF-8", warn = FALSE))
  if (length(genes) < 1) stop("诊断基因文件中没有有效基因。")
  genes
}

order_binary_group <- function(group) {
  group <- factor(group)
  if (nlevels(group) != 2) stop("诊断分析要求二分类；当前分组：", paste(levels(group), collapse = ", "))
  lev <- levels(group)
  control_hit <- grepl("^(con|control|normal|healthy)$", lev, ignore.case = TRUE)
  if (sum(control_hit) == 1) lev <- c(lev[control_hit], lev[!control_hit])
  factor(as.character(group), levels = lev)
}

make_stratified_folds <- function(y, k = 5, seed = 1L) {
  y <- factor(y)
  k <- min(k, min(table(y)))
  if (k < 3) stop("每组样本数太少，无法进行至少3折分层交叉验证。")
  set.seed(seed)
  foldid <- integer(length(y))
  for (lev in levels(y)) {
    idx <- which(y == lev)
    foldid[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  foldid
}

numeric_class_weights <- function(y_num) {
  tab <- table(y_num)
  unname(setNames(length(y_num) / (2 * as.numeric(tab)), names(tab))[as.character(y_num)])
}

# 岭回归(alpha=0)解决普通Logistic在小样本中的完全分离问题。
# 所有标准化和lambda选择均在训练折内完成，测试折只用于预测。
repeated_cv_ridge <- function(x, y_num, folds = 5L, repeats = 20L, seed = 1L) {
  x <- as.matrix(x)
  # glmnet要求至少两列；单基因模型加入恒为0的占位列，不影响预测或系数。
  if (ncol(x) == 1L) {
    x <- cbind(x, .glmnet_placeholder = 0)
  }
  pred_matrix <- matrix(NA_real_, nrow = nrow(x), ncol = repeats,
                        dimnames = list(rownames(x), paste0("Repeat", seq_len(repeats))))
  fold_records <- vector("list", repeats)
  for (r in seq_len(repeats)) {
    outer_id <- make_stratified_folds(factor(y_num, levels = c(0, 1)), folds, seed + r)
    fold_records[[r]] <- data.frame(Sample = rownames(x), Repeat = r, Fold = outer_id)
    for (fold in sort(unique(outer_id))) {
      test_idx <- which(outer_id == fold)
      train_idx <- which(outer_id != fold)
      y_train <- y_num[train_idx]
      inner_k <- min(3L, min(table(y_train)))
      inner_id <- make_stratified_folds(factor(y_train, levels = c(0, 1)), inner_k,
                                        seed + r * 100L + fold)
      train_weight <- numeric_class_weights(y_train)
      cv_fit <- cv.glmnet(
        x[train_idx, , drop = FALSE], y_train,
        family = "binomial", alpha = 0, standardize = TRUE,
        foldid = inner_id, nfolds = inner_k, weights = train_weight,
        type.measure = "deviance"
      )
      pred_matrix[test_idx, r] <- as.numeric(predict(
        cv_fit, newx = x[test_idx, , drop = FALSE], s = "lambda.1se", type = "response"
      ))
    }
  }
  list(prediction = rowMeans(pred_matrix, na.rm = TRUE),
       prediction_by_repeat = pred_matrix,
       folds = do.call(rbind, fold_records))
}

metric_at_half <- function(y, probability) {
  pred <- as.integer(probability >= 0.5)
  sensitivity <- mean(pred[y == 1] == 1)
  specificity <- mean(pred[y == 0] == 0)
  accuracy <- mean(pred == y)
  c(Accuracy = accuracy, Sensitivity = sensitivity, Specificity = specificity,
    BalancedAccuracy = mean(c(sensitivity, specificity)))
}

net_benefit <- function(y, probability, thresholds) {
  n <- length(y)
  vapply(thresholds, function(pt) {
    pred <- probability >= pt
    tp <- sum(pred & y == 1)
    fp <- sum(pred & y == 0)
    tp / n - fp / n * pt / (1 - pt)
  }, numeric(1))
}

cat("==== 稳健诊断模型分析启动 ====\n")
if (!file.exists(input_expr)) stop("未找到表达矩阵：", input_expr)
feature_genes <- read_gene_list(gene_file)
dat_expr <- read.csv(input_expr, header = TRUE, row.names = 1, check.names = FALSE,
                     fileEncoding = "UTF-8-BOM")
rownames(dat_expr) <- gsub("-", "_", rownames(dat_expr))
found_genes <- intersect(feature_genes, rownames(dat_expr))
missing_genes <- setdiff(feature_genes, rownames(dat_expr))
if (length(found_genes) < 2) stop("表达矩阵中匹配到的诊断基因少于2个。")
if (length(missing_genes) > 0) warning("未匹配基因：", paste(missing_genes, collapse = ", "))

x <- t(as.matrix(dat_expr[found_genes, , drop = FALSE]))
storage.mode(x) <- "double"
if (anyNA(x)) stop("诊断基因表达矩阵中存在NA。")
group <- order_binary_group(sub("^.*_", "", rownames(x)))
control_label <- levels(group)[1]
disease_label <- levels(group)[2]
y <- as.integer(group == disease_label)
outer_folds <- min(outer_folds, min(table(group)))
model_vars <- colnames(x)

cat(sprintf("样本=%d；%s=%d；%s=%d；基因=%s\n",
            nrow(x), control_label, sum(y == 0), disease_label, sum(y == 1),
            paste(model_vars, collapse = ", ")))

# 保存进入模型的数据，GroupType保持0/1，便于下游复核。
filtered_output <- data.frame(x, GroupType = y, OriginalGroup = as.character(group), check.names = FALSE)
write.csv(filtered_output, file.path(result_dir, "Filtered_Expression_Matrix.csv"), quote = FALSE)

cat("运行20次重复、5折分层交叉验证的岭Logistic联合模型...\n")
combined_cv <- repeated_cv_ridge(x, y, outer_folds, outer_repeats, base_seed)
combined_prob <- combined_cv$prediction
write.csv(combined_cv$folds, file.path(result_dir, "Diagnostic_Model_Stratified_Folds.csv"), row.names = FALSE)

# 单基因折外预测用于DCA；原始表达值仍用于传统单基因ROC。
gene_cv_prob <- matrix(NA_real_, nrow = nrow(x), ncol = ncol(x),
                       dimnames = list(rownames(x), model_vars))
for (i in seq_along(model_vars)) {
  cat(sprintf("运行单基因折外预测：%s (%d/%d)\n", model_vars[i], i, length(model_vars)))
  gene_cv_prob[, i] <- repeated_cv_ridge(
    x[, i, drop = FALSE], y, outer_folds, outer_repeats, base_seed + i * 1000L
  )$prediction
}

oof_df <- data.frame(Sample = rownames(x), Group = as.character(group), Outcome = y,
                     Combined_Model = combined_prob, gene_cv_prob, check.names = FALSE)
write.csv(oof_df, file.path(result_dir, "OOF_Predictions.csv"), row.names = FALSE)

# 在全数据上拟合最终岭模型，仅用于系数和列线图；性能评价一律使用上面的折外预测。
full_foldid <- make_stratified_folds(group, outer_folds, base_seed)
full_fit <- cv.glmnet(x, y, family = "binomial", alpha = 0, standardize = TRUE,
                      foldid = full_foldid, weights = numeric_class_weights(y),
                      type.measure = "deviance")
coef_mat <- as.matrix(coef(full_fit, s = "lambda.1se"))
coef_df <- data.frame(Term = rownames(coef_mat), Coefficient = as.numeric(coef_mat[, 1]),
                      Lambda = full_fit$lambda.1se, row.names = NULL)
write.csv(coef_df, file.path(result_dir, "Penalized_Model_Coefficients.csv"), row.names = FALSE)

# 兼容原输出名，但明确不为惩罚模型报告无效的普通回归P值。
coef_compat <- transform(coef_df, P_Value = NA_real_,
                         Note = "Ridge-penalized coefficient; ordinary Wald P-value not reported")
write.csv(coef_compat, file.path(result_dir, "Model_Significant_Coefficients.csv"), row.names = FALSE)

# 使用惩罚Logistic模型绘制列线图，解决完全分离造成的无限系数。
nomogram_ok <- FALSE
if (requireNamespace("rms", quietly = TRUE)) {
  nomo_data <- data.frame(x, GroupType = factor(group, levels = c(control_label, disease_label)),
                          check.names = FALSE)
  ddinfo <- rms::datadist(nomo_data)
  old_dd <- getOption("datadist")
  options(datadist = "ddinfo")
  assign("ddinfo", ddinfo, envir = .GlobalEnv)
  nomo_formula <- as.formula(paste("GroupType ~", paste(model_vars, collapse = " + ")))
  nomogram_ok <- tryCatch({
    nomo_fit <- rms::lrm(nomo_formula, data = nomo_data, x = TRUE, y = TRUE, penalty = 1)
    nomo_obj <- rms::nomogram(nomo_fit, fun = plogis,
                              fun.at = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95),
                              lp = FALSE, funlabel = paste0("Probability of ", disease_label))
    pdf(file.path(result_dir, "Nomogram_Plot.pdf"), width = 11, height = 6)
    plot(nomo_obj, xfrac = 0.35)
    dev.off()
    TRUE
  }, error = function(e) {
    warning("rms列线图生成失败，改为输出惩罚系数图：", e$message)
    FALSE
  })
  options(datadist = old_dd)
}
if (!nomogram_ok) {
  plot_coef <- coef_df[coef_df$Term != "(Intercept)", , drop = FALSE]
  pdf(file.path(result_dir, "Nomogram_Plot.pdf"), width = 8, height = 5)
  barplot(plot_coef$Coefficient, names.arg = plot_coef$Term, horiz = TRUE, las = 1,
          col = ifelse(plot_coef$Coefficient >= 0, "#D55E00", "#0072B2"),
          xlab = "Ridge-penalized coefficient", main = "Penalized Diagnostic Model")
  abline(v = 0, col = "gray40")
  dev.off()
}

# 联合模型折外ROC与评价指标。
combined_roc <- roc(y, combined_prob, levels = c(0, 1), direction = "<", quiet = TRUE)
combined_ci <- as.numeric(ci.auc(combined_roc))
combined_metrics <- metric_at_half(y, combined_prob)
combined_result <- data.frame(
  Model = "Ridge logistic (repeated stratified OOF)",
  AUC = as.numeric(auc(combined_roc)), AUC_CI_Lower = combined_ci[1], AUC_CI_Upper = combined_ci[3],
  Threshold = 0.5, Accuracy = combined_metrics["Accuracy"],
  Sensitivity = combined_metrics["Sensitivity"], Specificity = combined_metrics["Specificity"],
  BalancedAccuracy = combined_metrics["BalancedAccuracy"],
  BrierScore = mean((combined_prob - y)^2), OuterFolds = outer_folds,
  Repeats = outer_repeats, Validation = "Internal repeated stratified cross-validation"
)
write.csv(combined_result, file.path(result_dir, "Combined_Model_CV_Metrics.csv"), row.names = FALSE)

pdf(file.path(result_dir, "Combined_Model_CV_ROC.pdf"), width = 6.5, height = 6.5)
plot(combined_roc, legacy.axes = TRUE, col = "#D55E00", lwd = 3,
     main = "Combined Model ROC (Out-of-fold Predictions)")
abline(a = 0, b = 1, lty = 2, col = "gray65")
legend("bottomright",
       legend = sprintf("AUC = %.3f (95%% CI %.3f-%.3f)", combined_result$AUC,
                        combined_result$AUC_CI_Lower, combined_result$AUC_CI_Upper),
       col = "#D55E00", lwd = 3, bty = "n")
dev.off()

# 单基因ROC：用于描述每个基因自身的判别能力。
gene_colors <- c("#0072B2", "#009E73", "#E69F00", "#CC79A7", "#56B4E9", "#D55E00")
gene_colors <- rep(gene_colors, length.out = length(model_vars))
gene_rocs <- vector("list", length(model_vars)); names(gene_rocs) <- model_vars
gene_auc_rows <- vector("list", length(model_vars))
for (i in seq_along(model_vars)) {
  g <- model_vars[i]
  r <- roc(y, x[, g], levels = c(0, 1), direction = "auto", quiet = TRUE)
  ci_value <- as.numeric(ci.auc(r))
  mean_control <- mean(x[y == 0, g])
  mean_disease <- mean(x[y == 1, g])
  gene_rocs[[i]] <- r
  gene_auc_rows[[i]] <- data.frame(
    Gene = g, AUC = as.numeric(auc(r)), AUC_CI_Lower = ci_value[1], AUC_CI_Upper = ci_value[3],
    Direction = ifelse(mean_disease > mean_control,
                       paste0(disease_label, " higher"), paste0(control_label, " higher")),
    Mean_Control = mean_control, Mean_Disease = mean_disease
  )
}
gene_auc_df <- do.call(rbind, gene_auc_rows)
write.csv(gene_auc_df, file.path(result_dir, "IndividualGenes_AUC.csv"), row.names = FALSE)

pdf(file.path(result_dir, "IndividualGenes_ROC_SCI.pdf"), width = 7.5, height = 7.5)
plot(gene_rocs[[1]], legacy.axes = TRUE, col = gene_colors[1], lwd = 3,
     main = "Individual Gene ROC Curves")
if (length(gene_rocs) > 1) {
  for (i in 2:length(gene_rocs)) plot(gene_rocs[[i]], add = TRUE, col = gene_colors[i], lwd = 3)
}
abline(a = 0, b = 1, lty = 2, col = "gray70")
legend("bottomright", legend = sprintf("%s  AUC=%.3f", gene_auc_df$Gene, gene_auc_df$AUC),
       col = gene_colors, lwd = 3, bty = "n", cex = 0.95)
dev.off()

# 折外校准曲线：不再用训练集拟合值评价训练集。
clip_prob <- pmin(pmax(combined_prob, 1e-6), 1 - 1e-6)
calibration_fit <- tryCatch(glm(y ~ qlogis(clip_prob), family = binomial()), error = function(e) NULL)
calibration_intercept <- if (is.null(calibration_fit)) NA_real_ else unname(coef(calibration_fit)[1])
calibration_slope <- if (is.null(calibration_fit)) NA_real_ else unname(coef(calibration_fit)[2])
breaks <- unique(quantile(combined_prob, probs = seq(0, 1, length.out = 6), na.rm = TRUE))
if (length(breaks) < 3) breaks <- seq(0, 1, length.out = 4)
cal_bin <- cut(combined_prob, breaks = breaks, include.lowest = TRUE, labels = FALSE)
cal_data <- aggregate(cbind(Predicted = combined_prob, Observed = y), list(Bin = cal_bin), mean)
write.csv(cal_data, file.path(result_dir, "Calibration_Data.csv"), row.names = FALSE)
write.csv(data.frame(CalibrationIntercept = calibration_intercept,
                     CalibrationSlope = calibration_slope,
                     BrierScore = mean((combined_prob - y)^2)),
          file.path(result_dir, "Calibration_Metrics.csv"), row.names = FALSE)

pdf(file.path(result_dir, "Calibration_Curve.pdf"), width = 6, height = 6)
plot(cal_data$Predicted, cal_data$Observed, type = "b", pch = 19, lwd = 2,
     xlim = c(0, 1), ylim = c(0, 1), xlab = "Mean predicted probability",
     ylab = "Observed disease proportion", main = "OOF Calibration Curve")
abline(a = 0, b = 1, lty = 2, col = "gray55")
legend("topleft", legend = c("Repeated-CV OOF", "Ideal"),
       col = c("black", "gray55"), lty = c(1, 2), pch = c(19, NA), bty = "n")
dev.off()

# DCA完全基于折外预测，避免训练集拟合概率造成虚假的净获益=1。
thresholds <- seq(0.01, 0.99, by = 0.01)
prevalence <- mean(y)
dca_df <- data.frame(
  Threshold = thresholds,
  Diagnostic_Model = net_benefit(y, combined_prob, thresholds),
  Treat_All = prevalence - (1 - prevalence) * thresholds / (1 - thresholds),
  Treat_None = 0
)
for (g in model_vars) dca_df[[g]] <- net_benefit(y, gene_cv_prob[, g], thresholds)
write.csv(dca_df, file.path(result_dir, "Decision_Curve_Data.csv"), row.names = FALSE)

pdf(file.path(result_dir, "DCA.pdf"), width = 8, height = 7)
ylim_dca <- range(c(dca_df$Diagnostic_Model, dca_df$Treat_All, 0), finite = TRUE)
ylim_dca[1] <- max(ylim_dca[1], -0.25)
plot(thresholds, dca_df$Diagnostic_Model, type = "l", lwd = 3, col = "#D55E00",
     xlab = "Threshold Probability", ylab = "Net Benefit",
     main = "Decision Curve Analysis (OOF Predictions)", ylim = ylim_dca)
lines(thresholds, dca_df$Treat_All, lwd = 2, col = "black")
lines(thresholds, dca_df$Treat_None, lwd = 2, lty = 2, col = "gray45")
for (i in seq_along(model_vars)) lines(thresholds, dca_df[[model_vars[i]]], col = gene_colors[i], lwd = 1.5)
legend("topright", legend = c("Diagnostic Model", model_vars, "Treat All", "Treat None"),
       col = c("#D55E00", gene_colors, "black", "gray45"),
       lty = c(rep(1, length(model_vars) + 2), 2),
       lwd = c(3, rep(1.5, length(model_vars)), 2, 2), bty = "n", cex = 0.85)
dev.off()

# 小样本且不平衡，采用Wilcoxon秩和检验；箱线图按基因分面并使用独立y轴。
p_values <- vapply(model_vars, function(g) wilcox.test(x[y == 0, g], x[y == 1, g], exact = FALSE)$p.value,
                   numeric(1))
p_labels <- ifelse(p_values < 0.001, "***",
                   ifelse(p_values < 0.01, "**", ifelse(p_values < 0.05, "*", "ns")))
p_value_df <- data.frame(Gene = model_vars, P_Value = p_values, Significance = p_labels,
                         Test = "Wilcoxon rank-sum")
write.csv(p_value_df, file.path(result_dir, "Diagnostic_Genes_PValues.csv"), row.names = FALSE)

box_df <- do.call(rbind, lapply(model_vars, function(g) {
  data.frame(Gene = g, Expression = x[, g], Group = factor(group, levels = c(control_label, disease_label)))
}))
annotation_df <- do.call(rbind, lapply(seq_along(model_vars), function(i) {
  values <- x[, model_vars[i]]
  data.frame(Gene = model_vars[i], Group = disease_label,
             y = max(values) + max(diff(range(values)) * 0.10, 0.05), label = p_labels[i])
}))
p_box <- ggplot(box_df, aes(Group, Expression, fill = Group)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.78) +
  geom_jitter(width = 0.13, size = 1.5, alpha = 0.72) +
  geom_text(data = annotation_df, aes(x = Group, y = y, label = label),
            inherit.aes = FALSE, size = 5, fontface = "bold") +
  facet_wrap(~Gene, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = setNames(c("#B9D7EA", "#E98B8B"), c(control_label, disease_label))) +
  labs(x = NULL, y = "Expression level", title = "Diagnostic Gene Expression") +
  theme_classic(base_size = 12) +
  theme(legend.position = "top", strip.background = element_blank(),
        strip.text = element_text(face = "bold"), axis.text.x = element_text(angle = 35, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0.5))
ggsave(file.path(result_dir, "Diagnostic_Genes_Boxplot.pdf"), p_box, width = 10, height = 5.5)

analysis_notes <- c(
  "Model: ridge-penalized logistic regression (alpha=0).",
  sprintf("Validation: %d repeats of %d-fold stratified cross-validation.", outer_repeats, outer_folds),
  "ROC, calibration and DCA for the combined model use out-of-fold predictions.",
  "Class imbalance is handled with inverse-frequency weights.",
  "This remains internal validation; an independent external cohort is required before diagnostic claims.",
  "Upstream feature selection used the same discovery cohort, so some selection bias may remain."
)
writeLines(analysis_notes, file.path(result_dir, "Analysis_Notes.txt"), useBytes = TRUE)

cat("==== 分析完成 ====\n")
cat(sprintf("联合模型折外AUC：%.3f（95%%CI %.3f-%.3f）\n",
            combined_result$AUC, combined_result$AUC_CI_Lower, combined_result$AUC_CI_Upper))
cat(sprintf("平衡准确率：%.3f；Brier分数：%.3f\n",
            combined_result$BalancedAccuracy, combined_result$BrierScore))
cat("结果目录：", result_dir, "\n", sep = "")
