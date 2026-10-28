# -*- coding: utf-8 -*-
options(encoding = "UTF-8", stringsAsFactors = FALSE)
try(Sys.setlocale(category = "LC_CTYPE", locale = "Chinese_China.65001"), silent = TRUE)

suppressPackageStartupMessages({
  library(glmnet)
  library(reshape2)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
})

# 本脚本直接对应当前项目目录。若以后移动项目，只需修改这一行。
workDir <- Sys.getenv("IPF_BULK_WORKDIR", unset = getwd())
workDir <- normalizePath(workDir, winslash = "/", mustWork = TRUE)
setwd(workDir)

inputFile <- file.path(workDir, "Sample Type Matrix.csv")
geneListFile <- file.path(workDir, "Intersection_DEG_ferro-aging.csv")
outputDir <- file.path(workDir, "output")
if (!dir.exists(outputDir)) dir.create(outputDir, recursive = TRUE)

strip_bom <- function(x) sub("^\ufeff", "", x)
clean_gene_vector <- function(x) {
  x <- strip_bom(trimws(as.character(x)))
  x <- gsub("^['\"]|['\"]$", "", x)
  unique(x[!is.na(x) & x != ""])
}

read_gene_list <- function(file_path) {
  if (!file.exists(file_path)) stop("未找到待研究基因文件：", file_path)
  gene_df <- read.csv(file_path, header = TRUE, check.names = FALSE,
                      stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  names(gene_df) <- strip_bom(names(gene_df))
  gene_col <- which(tolower(gsub("[^A-Za-z0-9]+", "", names(gene_df))) %in%
                      c("gene", "genes", "symbol", "genesymbol"))
  if (length(gene_col) == 0) gene_col <- 1L
  genes <- clean_gene_vector(gene_df[[gene_col[1]]])
  if (length(genes) == 0) stop("待研究基因文件中没有有效基因：", file_path)
  genes
}

make_stratified_folds <- function(y, k = 5, seed = 123456) {
  y <- factor(y)
  k <- min(k, min(table(y)))
  if (k < 3) stop("每组样本数太少，无法进行至少3折的分层交叉验证。")
  set.seed(seed)
  foldid <- integer(length(y))
  for (lev in levels(y)) {
    idx <- which(y == lev)
    foldid[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  list(foldid = foldid, k = k)
}

order_binary_group <- function(group) {
  group <- factor(group)
  if (nlevels(group) != 2) {
    stop("样本分组必须是二分类；当前识别到：", paste(levels(group), collapse = ", "))
  }
  lev <- levels(group)
  control_hit <- grepl("^(con|control|normal|healthy)$", lev, ignore.case = TRUE)
  if (sum(control_hit) == 1) lev <- c(lev[control_hit], lev[!control_hit])
  factor(as.character(group), levels = lev)
}

message("读取表达矩阵：", inputFile)
if (!file.exists(inputFile)) stop("未找到表达矩阵文件：", inputFile)
expr_df <- read.csv(inputFile, header = TRUE, check.names = FALSE,
                    stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
if (ncol(expr_df) < 3) stop("表达矩阵列数不足：第一列应为基因名，其余列为样本。")
names(expr_df) <- strip_bom(names(expr_df))

gene_symbols <- strip_bom(trimws(as.character(expr_df[[1]])))
keep <- !is.na(gene_symbols) & gene_symbols != ""
gene_symbols <- gene_symbols[keep]
expr_df <- expr_df[keep, -1, drop = FALSE]
expr_df[] <- lapply(expr_df, function(z) suppressWarnings(as.numeric(z)))
dedup <- !duplicated(gene_symbols)
expr_df <- expr_df[dedup, , drop = FALSE]
gene_symbols <- gene_symbols[dedup]
rownames(expr_df) <- gene_symbols
expr_df <- expr_df[rowSums(is.na(expr_df)) < ncol(expr_df), , drop = FALSE]
expr_df[is.na(expr_df)] <- 0

study_genes <- read_gene_list(geneListFile)
available_genes <- intersect(study_genes, rownames(expr_df))
missing_genes <- setdiff(study_genes, rownames(expr_df))
if (length(available_genes) < 2) stop("匹配到的待研究基因少于2个，无法进行LASSO。")
write.csv(data.frame(Gene = available_genes),
          file.path(outputDir, "Available_Study_Genes.csv"), row.names = FALSE)
write.csv(data.frame(Gene = missing_genes),
          file.path(outputDir, "Missing_Study_Genes.csv"), row.names = FALSE)

x_matrix <- t(as.matrix(expr_df[available_genes, , drop = FALSE]))
sample_groups <- order_binary_group(sub("^.*_", "", rownames(x_matrix)))
y_numeric <- as.integer(sample_groups == levels(sample_groups)[2])
group_tab <- table(sample_groups)
message("样本分组：", paste(names(group_tab), as.integer(group_tab), sep = "=", collapse = "; "))

# 5折分层交叉验证：8例对照时，每折都能保留对照样本。
fold_obj <- make_stratified_folds(sample_groups, k = 5, seed = 123456)
foldid <- fold_obj$foldid
nfolds <- fold_obj$k
fold_table <- as.data.frame.matrix(table(Fold = foldid, Group = sample_groups))
fold_table$Fold <- rownames(fold_table)
rownames(fold_table) <- NULL
write.csv(fold_table, file.path(outputDir, "LASSO_Fold_Distribution.csv"), row.names = FALSE)

# 逆频率权重，降低40:8类别不平衡对模型的影响。
class_weight <- setNames(length(sample_groups) / (2 * as.numeric(group_tab)), names(group_tab))
obs_weight <- unname(class_weight[as.character(sample_groups)])

set.seed(123456)
lasso_model <- glmnet(x_matrix, y_numeric, family = "binomial", alpha = 1,
                      standardize = TRUE, weights = obs_weight)
cv_model <- cv.glmnet(
  x_matrix, y_numeric,
  family = "binomial", alpha = 1, type.measure = "deviance",
  foldid = foldid, nfolds = nfolds, standardize = TRUE,
  weights = obs_weight, keep = TRUE
)

# lambda.min用于候选基因发现；lambda.1se结果同时保存，供更保守选择时参考。
extract_coef <- function(s_value) {
  co <- coef(cv_model, s = s_value)
  values <- as.numeric(co[-1, 1])
  genes <- rownames(co)[-1]
  keep <- values != 0
  out <- data.frame(gene = genes[keep], coefficient = values[keep],
                    abs_coefficient = abs(values[keep]), stringsAsFactors = FALSE)
  out[order(out$abs_coefficient, decreasing = TRUE), , drop = FALSE]
}
selected_min <- extract_coef("lambda.min")
selected_1se <- extract_coef("lambda.1se")
write.csv(selected_min, file.path(outputDir, "LASSO.gene.csv"), row.names = FALSE, quote = FALSE)
write.csv(selected_1se, file.path(outputDir, "LASSO.gene.lambda_1se.csv"), row.names = FALSE, quote = FALSE)
if (nrow(selected_min) > 0) {
  write.csv(expr_df[selected_min$gene, , drop = FALSE],
            file.path(outputDir, "LASSO.gene.exp.csv"), quote = FALSE)
}

gene_count <- vapply(cv_model$lambda, function(lam) {
  sum(as.numeric(coef(lasso_model, s = lam)[-1, 1]) != 0)
}, integer(1))
cv_data <- data.frame(lambda = cv_model$lambda, log_lambda = log(cv_model$lambda),
                      cvm = cv_model$cvm, cvup = cv_model$cvup,
                      cvlo = cv_model$cvlo, gene_count = gene_count)
write.csv(cv_data, file.path(outputDir, "cv_results_table_enhanced.csv"), row.names = FALSE)

all_coef_df <- as.data.frame(as.matrix(lasso_model$beta))
all_coef_df$gene <- rownames(all_coef_df)
all_coef_df <- all_coef_df[, c("gene", setdiff(names(all_coef_df), "gene")), drop = FALSE]
write.csv(all_coef_df, file.path(outputDir, "all_coefficients_table_enhanced.csv"), row.names = FALSE)

lambda_min_count <- nrow(selected_min)
lambda_1se_count <- nrow(selected_1se)
p_count <- ggplot(cv_data, aes(gene_count, cvm)) +
  geom_errorbar(aes(ymin = cvlo, ymax = cvup), width = 0.3, color = "#7F7F7F") +
  geom_point(size = 2.5, color = "#0072B2") +
  geom_vline(xintercept = c(lambda_min_count, lambda_1se_count), linetype = "dashed", color = "red") +
  labs(x = "Number of Selected Genes", y = "Mean Cross-Validated Deviance",
       title = "Stratified 5-fold CV for LASSO") +
  theme_bw(base_size = 13) + theme(plot.title = element_text(face = "bold", hjust = 0.5))
ggsave(file.path(outputDir, "CV_curve_plot_A_gene_count.pdf"), p_count, width = 6, height = 5.5)

p_lambda <- ggplot(cv_data, aes(log_lambda, cvm)) +
  geom_errorbar(aes(ymin = cvlo, ymax = cvup), width = 0.05, color = "#7F7F7F") +
  geom_point(size = 2.5, color = "#0072B2") +
  geom_vline(xintercept = log(c(cv_model$lambda.min, cv_model$lambda.1se)),
             linetype = "dashed", color = "red") +
  labs(x = "Log(lambda)", y = "Mean Cross-Validated Deviance",
       title = "Stratified 5-fold CV for LASSO") +
  theme_bw(base_size = 13) + theme(plot.title = element_text(face = "bold", hjust = 0.5))
ggsave(file.path(outputDir, "CV_curve_plot_A.pdf"), p_lambda, width = 6, height = 5.5)

lasso_long <- reshape2::melt(
  data.frame(log_lambda = log(lasso_model$lambda), t(as.matrix(lasso_model$beta)), check.names = FALSE),
  id.vars = "log_lambda"
)
representative_df <- do.call(rbind, lapply(split(lasso_long, lasso_long$variable), function(z) {
  z[which.max(abs(z$value)), , drop = FALSE]
}))
n_vars <- length(unique(lasso_long$variable))
palette <- colorRampPalette(brewer.pal(11, "Spectral"))(max(n_vars, 3))
p_path <- ggplot(lasso_long, aes(log_lambda, value, group = variable, color = variable)) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  geom_text_repel(data = representative_df, aes(label = variable), size = 3,
                  max.overlaps = Inf, show.legend = FALSE) +
  geom_vline(xintercept = log(c(cv_model$lambda.min, cv_model$lambda.1se)),
             linetype = "dashed", color = "#D55E00") +
  scale_color_manual(values = palette[seq_len(n_vars)]) +
  labs(x = "Log(lambda)", y = "Coefficient", title = "LASSO Coefficient Path") +
  theme_bw(base_size = 13) + theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5))
ggsave(file.path(outputDir, "LASSO_coefficient_path_plot_B.pdf"), p_path, width = 8, height = 6)

summary_df <- data.frame(
  Item = c("input_matrix", "study_gene_file", "study_genes", "matched_genes", "samples",
           "groups", "validation", "class_weighting", "lambda.min", "lambda.1se",
           "selected_genes_lambda_min", "selected_genes_lambda_1se"),
  Value = c(inputFile, geneListFile, length(study_genes), length(available_genes), nrow(x_matrix),
            paste(names(group_tab), as.integer(group_tab), sep = "=", collapse = "; "),
            paste0(nfolds, "-fold stratified CV"), "inverse-frequency",
            cv_model$lambda.min, cv_model$lambda.1se, nrow(selected_min), nrow(selected_1se))
)
write.csv(summary_df, file.path(outputDir, "LASSO_Run_Summary.csv"), row.names = FALSE)

cat("LASSO分析完成\n")
cat("分层交叉验证折数：", nfolds, "\n", sep = "")
cat("lambda.min筛选基因数：", nrow(selected_min), "\n", sep = "")
cat("lambda.1se筛选基因数：", nrow(selected_1se), "\n", sep = "")
cat("输出目录：", outputDir, "\n", sep = "")
