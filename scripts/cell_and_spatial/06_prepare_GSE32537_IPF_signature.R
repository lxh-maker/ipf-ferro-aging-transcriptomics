suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
root <- normalizePath(file.path(dirname(script_file), "..", ".."), mustWork = FALSE)
data_dir <- file.path(root, "00_data", "GSE32537")
out_dir <- file.path(root, "05_drug_screen", "GSE32537_signature")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

matrix_file <- file.path(data_dir, "GSE32537_series_matrix.txt.gz")
annotation_file <- file.path(data_dir, "GPL6244.annot.gz")

unquote <- function(x) gsub('^"|"$', "", x)
split_geo_line <- function(line) unquote(strsplit(line, "\t", fixed = TRUE)[[1]][-1])

message("Reading GSE32537 series matrix")
con <- gzfile(matrix_file, open = "rt")
matrix_lines <- readLines(con, warn = FALSE)
close(con)
table_begin <- match("!series_matrix_table_begin", matrix_lines)
table_end <- match("!series_matrix_table_end", matrix_lines)
meta_lines <- matrix_lines[seq_len(table_begin - 1L)]
sample_ids <- split_geo_line(meta_lines[grepl("^!Sample_geo_accession", meta_lines)][1])
diagnosis_line <- meta_lines[
  grepl("^!Sample_characteristics_ch1", meta_lines) &
    grepl("final diagnosis:", meta_lines, ignore.case = TRUE)
][1]
diagnosis <- sub("^final diagnosis:\\s*", "", split_geo_line(diagnosis_line), ignore.case = TRUE)

tc <- textConnection(paste(matrix_lines[(table_begin + 1L):(table_end - 1L)], collapse = "\n"))
expr_tab <- read.delim(tc, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, quote = '"')
close(tc)
rm(matrix_lines); invisible(gc())
rownames(expr_tab) <- as.character(expr_tab[[1]])
expr <- as.matrix(expr_tab[, -1, drop = FALSE])
storage.mode(expr) <- "double"
colnames(expr) <- sample_ids
rm(expr_tab); invisible(gc())

message("Reading GPL6244 annotation")
con <- gzfile(annotation_file, open = "rt")
annot_lines <- readLines(con, warn = FALSE)
close(con)
annot_begin <- match("!platform_table_begin", annot_lines)
annot_end <- match("!platform_table_end", annot_lines)
tc <- textConnection(paste(annot_lines[(annot_begin + 1L):(annot_end - 1L)], collapse = "\n"))
annotation <- read.delim(tc, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, quote = '"')
close(tc)
rm(annot_lines); invisible(gc())
stopifnot(all(c("ID", "Gene symbol") %in% colnames(annotation)))

# Select one probe per HGNC symbol using the largest across-cohort IQR. This
# unsupervised choice does not use the IPF/control labels.
probe_iqr <- apply(expr, 1, IQR, na.rm = TRUE)
annotation <- annotation[annotation$ID %in% rownames(expr), , drop = FALSE]
symbol_lists <- strsplit(as.character(annotation[["Gene symbol"]]), " /// ", fixed = TRUE)
map <- do.call(rbind, lapply(seq_along(symbol_lists), function(i) {
  z <- unique(trimws(symbol_lists[[i]]))
  z <- z[nzchar(z) & z != "---"]
  if (!length(z)) return(NULL)
  data.frame(Gene = z, Probe = as.character(annotation$ID[i]), stringsAsFactors = FALSE)
}))
map$IQR <- unname(probe_iqr[map$Probe])
map <- map[is.finite(map$IQR) & grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", map$Gene), , drop = FALSE]
map <- map[order(map$Gene, -map$IQR, map$Probe), , drop = FALSE]
chosen <- map[!duplicated(map$Gene), , drop = FALSE]
gene_expr <- expr[chosen$Probe, , drop = FALSE]
rownames(gene_expr) <- chosen$Gene
rm(expr); invisible(gc())

keep <- diagnosis %in% c("control", "IPF/UIP")
gene_expr <- gene_expr[, keep, drop = FALSE]
group <- factor(diagnosis[keep], levels = c("control", "IPF/UIP"))
sample_meta <- data.frame(Sample = colnames(gene_expr), Group = as.character(group), stringsAsFactors = FALSE)
write.csv(sample_meta, file.path(out_dir, "GSE32537_sample_metadata.csv"), row.names = FALSE)
write.csv(chosen, file.path(out_dir, "GSE32537_selected_probe_per_gene.csv"), row.names = FALSE)
saveRDS(gene_expr, file.path(out_dir, "GSE32537_gene_expression.rds"), compress = "xz")

design <- model.matrix(~0 + group)
colnames(design) <- c("Control", "IPF")
contrast <- makeContrasts(IPF - Control, levels = design)
fit <- eBayes(contrasts.fit(lmFit(gene_expr, design), contrast), robust = TRUE)
de <- topTable(fit, number = Inf, sort.by = "P")
de$Gene <- rownames(de)
de <- de[, c("Gene", setdiff(colnames(de), "Gene")), drop = FALSE]
de$Regulation <- ifelse(de$adj.P.Val < 0.05 & de$logFC > 0.5, "Up",
                        ifelse(de$adj.P.Val < 0.05 & de$logFC < -0.5, "Down", "Not_significant"))
write.csv(de, file.path(out_dir, "GSE32537_IPF_vs_Control_limma_all.csv"), row.names = FALSE)
write.csv(de[de$Regulation != "Not_significant", ],
          file.path(out_dir, "GSE32537_IPF_vs_Control_DEG_FDR05_logFC05.csv"), row.names = FALSE)

control_idx <- group == "control"
ipf_idx <- group == "IPF/UIP"
control_mean <- rowMeans(gene_expr[, control_idx, drop = FALSE], na.rm = TRUE)
ipf_mean <- rowMeans(gene_expr[, ipf_idx, drop = FALSE], na.rm = TRUE)
control_var <- apply(gene_expr[, control_idx, drop = FALSE], 1, var, na.rm = TRUE)
ipf_var <- apply(gene_expr[, ipf_idx, drop = FALSE], 1, var, na.rm = TRUE)
denominator <- sqrt(control_var + ipf_var)
denominator[!is.finite(denominator) | denominator == 0] <- 1e-10
vscore <- (ipf_mean - control_mean) / denominator
vtab <- data.frame(Gene = names(vscore), vscore = as.numeric(vscore),
                   logFC = de$logFC[match(names(vscore), de$Gene)],
                   adj.P.Val = de$adj.P.Val[match(names(vscore), de$Gene)])
vtab <- vtab[order(-abs(vtab$vscore)), ]
write.csv(vtab, file.path(out_dir, "GSE32537_IPF_vs_Control_vscores.csv"), row.names = FALSE)

plot_de <- de
plot_de$minus_log10_fdr <- -log10(pmax(plot_de$adj.P.Val, .Machine$double.xmin))
p <- ggplot(plot_de, aes(logFC, minus_log10_fdr, color = Regulation)) +
  geom_point(alpha = 0.55, size = 1.15) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = 2, color = "grey45") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey45") +
  scale_color_manual(values = c(Down = "#2563EB", Not_significant = "grey75", Up = "#DC2626")) +
  theme_bw(base_size = 12) +
  labs(title = "GSE32537: IPF/UIP versus control lung", x = "log2 fold change", y = "-log10(FDR)", color = NULL)
ggsave(file.path(out_dir, "GSE32537_IPF_vs_Control_volcano.pdf"), p, width = 7, height = 5.5)
ggsave(file.path(out_dir, "GSE32537_IPF_vs_Control_volcano.png"), p, width = 7, height = 5.5, dpi = 300)

summary_lines <- c(
  sprintf("Samples: %d control, %d IPF/UIP", sum(control_idx), sum(ipf_idx)),
  sprintf("Mapped genes: %d", nrow(gene_expr)),
  sprintf("DEGs (FDR<0.05, |logFC|>0.5): %d up, %d down",
          sum(de$Regulation == "Up"), sum(de$Regulation == "Down")),
  "Contrast: IPF/UIP - control",
  "Probe selection: maximum across-cohort IQR without outcome labels",
  "V-score: (mean_IPF - mean_control) / sqrt(var_IPF + var_control)"
)
writeLines(summary_lines, file.path(out_dir, "SIGNATURE_SUMMARY.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"), useBytes = TRUE)
cat(paste(summary_lines, collapse = "\n"), "\n")

