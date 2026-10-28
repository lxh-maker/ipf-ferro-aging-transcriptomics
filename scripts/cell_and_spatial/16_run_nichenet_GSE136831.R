#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(ggplot2)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
root <- normalizePath(file.path(dirname(script_file), "..", ".."), mustWork = FALSE)
data_dir <- file.path(root, "00_data", "NicheNet")
out_dir <- file.path(root, "06_mechanism_axis", "nichenet_GSE136831")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20261027)

read_pb <- function(category, type = "counts") {
  path <- file.path(out_dir, sprintf("pseudobulk_%s_%s.tsv.gz", type, category))
  x <- fread(path)
  genes <- x$gene
  mat <- as.matrix(x[, -"gene"])
  rownames(mat) <- genes
  storage.mode(mat) <- "numeric"
  mat
}

counts_at2 <- read_pb("ATII", "counts")
counts_sec <- read_pb("Secretory_Club_Goblet", "counts")
counts_mye <- read_pb("Myeloid", "counts")
detect_sec <- read_pb("Secretory_Club_Goblet", "detected")
detect_mye <- read_pb("Myeloid", "detected")

summary <- fread(file.path(out_dir, "donor_category_summary.csv"))
wide_cells <- dcast(summary, Subject_Identity ~ category, value.var = "n_cells", fill = 0)
paired_donors <- wide_cells[ATII >= 5 & Secretory_Club_Goblet >= 5, Subject_Identity]
if (length(paired_donors) != 15) stop("Expected 15 paired IPF donors, found ", length(paired_donors))

# Paired full-transcriptome pseudobulk differential expression.
pair_counts <- cbind(counts_at2[, paired_donors, drop = FALSE], counts_sec[, paired_donors, drop = FALSE])
colnames(pair_counts) <- c(paste0(paired_donors, "_ATII"), paste0(paired_donors, "_Secretory"))
state <- factor(rep(c("ATII", "Secretory"), each = length(paired_donors)), levels = c("ATII", "Secretory"))
donor <- factor(rep(paired_donors, 2))
design <- model.matrix(~ donor + state)

dge <- DGEList(counts = pair_counts)
keep <- filterByExpr(dge, design = design, min.count = 5)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
v <- voom(dge, design, plot = FALSE)
fit <- eBayes(lmFit(v, design), trend = FALSE, robust = TRUE)
de <- topTable(fit, coef = "stateSecretory", number = Inf, sort.by = "P")
de$gene <- rownames(de)
de <- de[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
fwrite(de, file.path(out_dir, "paired_secretory_vs_ATII_full_transcriptome_DE.csv"))

ligand_target_matrix <- readRDS(file.path(data_dir, "ligand_target_matrix_nsga2r_final.rds"))
lr_network <- readRDS(file.path(data_dir, "lr_network_human_21122021.rds"))

# Receiver background: genes retained by edgeR and present in the NicheNet model.
background <- intersect(rownames(dge), rownames(ligand_target_matrix))
de_bg <- de[match(background, de$gene), ]
de_bg <- de_bg[!is.na(de_bg$gene), ]

up_sig <- de_bg$gene[de_bg$adj.P.Val < 0.05 & de_bg$logFC > 1]
down_sig <- de_bg$gene[de_bg$adj.P.Val < 0.05 & de_bg$logFC < -1]
up_ranked <- de_bg[order(-de_bg$t), "gene"]
down_ranked <- de_bg[order(de_bg$t), "gene"]

geneset_up_250 <- head(intersect(up_ranked, up_sig), 250)
geneset_up_500 <- head(intersect(up_ranked, up_sig), 500)
geneset_up_all <- up_sig
geneset_down_500 <- head(intersect(down_ranked, down_sig), 500)

if (length(geneset_up_500) < 100) stop("Too few upregulated receiver genes for NicheNet")

# Expression filters: a candidate ligand must be detected in >=5% of IPF myeloid
# cells and have at least one LR-network receptor detected in >=5% of secretory cells.
cell_totals <- dcast(summary, Subject_Identity ~ category, value.var = "n_cells", fill = 0)
total_mye_cells <- sum(cell_totals$Myeloid)
total_sec_cells <- sum(cell_totals$Secretory_Club_Goblet)
ligand_detect_fraction <- rowSums(detect_mye) / total_mye_cells
ligand_detect_fraction_secretory <- rowSums(detect_sec) / total_sec_cells
receptor_detect_fraction <- rowSums(detect_sec) / total_sec_cells

expressed_ligands <- names(ligand_detect_fraction)[ligand_detect_fraction >= 0.05]
expressed_receptors <- names(receptor_detect_fraction)[receptor_detect_fraction >= 0.05]
lr_expressed <- unique(lr_network[
  lr_network$from %in% expressed_ligands & lr_network$to %in% expressed_receptors,
  c("from", "to", "database", "source")
])
potential_ligands <- intersect(unique(lr_expressed$from), colnames(ligand_target_matrix))
if (length(potential_ligands) < 20) stop("Too few expression-supported candidate ligands")

aupr_trapezoid <- function(scores, response) {
  ord <- order(scores, decreasing = TRUE)
  y <- as.integer(response[ord])
  npos <- sum(y)
  if (npos == 0 || npos == length(y)) return(NA_real_)
  recall <- c(0, cumsum(y) / npos)
  precision <- c(1, cumsum(y) / seq_along(y))
  sum(diff(recall) * (head(precision, -1) + tail(precision, -1)) / 2)
}

score_geneset <- function(geneset, label) {
  bg <- intersect(background, rownames(ligand_target_matrix))
  response <- bg %in% geneset
  npos <- sum(response)
  nneg <- length(response) - npos
  rows <- lapply(potential_ligands, function(ligand) {
    scores <- ligand_target_matrix[bg, ligand]
    pearson <- suppressWarnings(cor(scores, as.numeric(response), method = "pearson"))
    ranks <- rank(scores, ties.method = "average")
    auroc <- (sum(ranks[response]) - npos * (npos + 1) / 2) / (npos * nneg)
    aupr <- aupr_trapezoid(scores, response)
    top250 <- names(sort(scores, decreasing = TRUE))[seq_len(min(250, length(scores)))]
    overlap <- sum(top250 %in% geneset)
    tab <- matrix(c(overlap, length(top250) - overlap,
                    npos - overlap, nneg - (length(top250) - overlap)), nrow = 2)
    fisher_p <- fisher.test(tab, alternative = "greater")$p.value
    data.frame(
      geneset = label, ligand = ligand, pearson = pearson, auroc = auroc,
      aupr = aupr, aupr_corrected = aupr - npos / length(bg),
      top250_overlap = overlap, fisher_p = fisher_p,
      ligand_detect_fraction_myeloid = ligand_detect_fraction[ligand],
      ligand_detect_fraction_secretory = ligand_detect_fraction_secretory[ligand],
      myeloid_vs_secretory_detection_log2ratio = log2((ligand_detect_fraction[ligand] + 0.01) /
                                                        (ligand_detect_fraction_secretory[ligand] + 0.01)),
      n_expressed_receptors = length(unique(lr_expressed$to[lr_expressed$from == ligand])),
      stringsAsFactors = FALSE
    )
  })
  result <- rbindlist(rows, fill = TRUE)
  result$fisher_FDR <- p.adjust(result$fisher_p, method = "BH")
  result$rank_pearson <- rank(-result$pearson, ties.method = "min")
  result$rank_aupr_corrected <- rank(-result$aupr_corrected, ties.method = "min")
  result[order(-pearson)]
}

activity <- rbindlist(list(
  score_geneset(geneset_up_250, "Secretory_up_top250"),
  score_geneset(geneset_up_500, "Secretory_up_top500_primary"),
  score_geneset(geneset_up_all, "Secretory_up_all_FDR05_logFC1"),
  score_geneset(geneset_down_500, "ATII_down_top500")
), fill = TRUE)
fwrite(activity, file.path(out_dir, "nichenet_ligand_activity_all_sensitivities.csv"))

primary <- activity[geneset == "Secretory_up_top500_primary"][order(-pearson)]
fwrite(primary, file.path(out_dir, "nichenet_ligand_activity_primary.csv"))
fwrite(lr_expressed, file.path(out_dir, "expression_supported_ligand_receptor_pairs.csv"))

gene_sets <- rbindlist(list(
  data.table(gene = geneset_up_250, geneset = "Secretory_up_top250"),
  data.table(gene = geneset_up_500, geneset = "Secretory_up_top500_primary"),
  data.table(gene = geneset_up_all, geneset = "Secretory_up_all_FDR05_logFC1"),
  data.table(gene = geneset_down_500, geneset = "ATII_down_top500")
))
fwrite(gene_sets, file.path(out_dir, "nichenet_receiver_gene_sets.csv"))

# Top active ligand-target links: top 250 NicheNet targets intersecting the primary gene set.
top_ligands <- head(primary$ligand, 20)
link_rows <- lapply(top_ligands, function(ligand) {
  scores <- ligand_target_matrix[background, ligand]
  target_order <- names(sort(scores, decreasing = TRUE))
  targets <- intersect(head(target_order, 250), geneset_up_500)
  if (length(targets) == 0) return(NULL)
  de_sub <- de[match(targets, de$gene), ]
  data.table(
    ligand = ligand, target = targets,
    regulatory_potential = ligand_target_matrix[targets, ligand],
    target_rank_for_ligand = match(targets, target_order),
    target_logFC = de_sub$logFC,
    target_FDR = de_sub$adj.P.Val
  )
})
links <- rbindlist(link_rows, fill = TRUE)
fwrite(links, file.path(out_dir, "top20_ligand_target_links_primary.csv"))

# Direct, auditable look at the proposed epithelial-axis genes, regardless of top-250 cutoff.
focus_ligands <- unique(c("SPP1", "S100A8", "S100A9", "TGFB1", "OSM", "LIF", "AREG", head(primary$ligand, 10)))
focus_ligands <- intersect(focus_ligands, colnames(ligand_target_matrix))
axis_genes <- intersect(c("LCN2", "ACSL4", "SCGB1A1", "SCGB3A1", "MUC5B", "SFTPC", "SFTPA1",
                          "KRT17", "KRT19", "HIF1A", "EGR1", "PTGS2"), rownames(ligand_target_matrix))
axis_rows <- lapply(focus_ligands, function(ligand) {
  scores_all <- ligand_target_matrix[, ligand]
  ordered <- names(sort(scores_all, decreasing = TRUE))
  data.table(
    ligand = ligand, target = axis_genes,
    regulatory_potential = ligand_target_matrix[axis_genes, ligand],
    global_target_rank = match(axis_genes, ordered),
    global_target_percentile = 1 - (match(axis_genes, ordered) - 1) / length(ordered),
    target_logFC = de$logFC[match(axis_genes, de$gene)],
    target_FDR = de$adj.P.Val[match(axis_genes, de$gene)]
  )
})
axis_links <- rbindlist(axis_rows)
fwrite(axis_links, file.path(out_dir, "focus_ligand_epithelial_axis_links.csv"))

# Summarize ranks of prespecified candidates across sensitivity definitions.
prespecified <- c("SPP1", "S100A8", "S100A9", "TGFB1", "TGFB2", "TGFB3", "OSM", "LIF", "AREG")
sensitivity <- activity[ligand %in% prespecified,
  .(geneset, ligand, pearson, aupr_corrected, top250_overlap, fisher_FDR,
    rank_pearson, rank_aupr_corrected, ligand_detect_fraction_myeloid, n_expressed_receptors)]
fwrite(sensitivity, file.path(out_dir, "prespecified_ligand_sensitivity_summary.csv"))

focused_ligands <- intersect(c("TGFB1", "S100A8", "S100A9", "SPP1", "OSM", "AREG"), primary$ligand)
focused_source <- primary[ligand %in% focused_ligands,
  .(ligand, pearson, aupr_corrected, rank_pearson, top250_overlap, fisher_FDR,
    ligand_detect_fraction_myeloid, ligand_detect_fraction_secretory,
    myeloid_vs_secretory_detection_log2ratio, n_expressed_receptors)]
fwrite(focused_source, file.path(out_dir, "focused_ligand_source_activity_evidence.csv"))
focused_lr <- as.data.table(lr_expressed)[from %in% focused_ligands]
focused_lr$receptor_detect_fraction_secretory <- receptor_detect_fraction[focused_lr$to]
focused_lr <- focused_lr[order(from, -receptor_detect_fraction_secretory)]
fwrite(focused_lr, file.path(out_dir, "focused_ligand_receptor_expression_evidence.csv"))

# Leave-one-donor-out robustness: re-estimate the paired receiver gene set 15
# times while keeping the expression-supported candidate-ligand universe fixed.
loo_rows <- lapply(paired_donors, function(dropped_donor) {
  use_donors <- setdiff(paired_donors, dropped_donor)
  loo_counts <- cbind(counts_at2[, use_donors, drop = FALSE], counts_sec[, use_donors, drop = FALSE])
  colnames(loo_counts) <- c(
    paste0(use_donors, "_ATII"),
    paste0(use_donors, "_Secretory")
  )
  loo_state <- factor(rep(c("ATII", "Secretory"), each = length(use_donors)),
                      levels = c("ATII", "Secretory"))
  loo_donor <- factor(rep(use_donors, 2))
  loo_design <- model.matrix(~ loo_donor + loo_state)
  loo_dge <- DGEList(counts = loo_counts)
  loo_keep <- filterByExpr(loo_dge, design = loo_design, min.count = 5)
  loo_dge <- loo_dge[loo_keep, , keep.lib.sizes = FALSE]
  loo_dge <- calcNormFactors(loo_dge)
  loo_v <- voom(loo_dge, loo_design, plot = FALSE)
  loo_fit <- eBayes(lmFit(loo_v, loo_design), robust = TRUE)
  loo_de <- topTable(loo_fit, coef = "loo_stateSecretory", number = Inf, sort.by = "P")
  loo_de$gene <- rownames(loo_de)
  loo_de <- loo_de[loo_de$gene %in% background, ]
  loo_sig <- loo_de$gene[loo_de$adj.P.Val < 0.05 & loo_de$logFC > 1]
  loo_ranked <- loo_de$gene[order(-loo_de$t)]
  loo_geneset <- head(intersect(loo_ranked, loo_sig), 500)
  if (length(loo_geneset) < 100) loo_geneset <- head(loo_ranked[loo_de$t[match(loo_ranked, loo_de$gene)] > 0], 500)
  scored <- score_geneset(loo_geneset, paste0("LOO_", dropped_donor))
  scored$dropped_donor <- dropped_donor
  scored$n_receiver_genes <- length(loo_geneset)
  scored
})
loo_activity <- rbindlist(loo_rows, fill = TRUE)
fwrite(loo_activity, file.path(out_dir, "nichenet_leave_one_donor_out_activity.csv"))
loo_summary <- loo_activity[, .(
  median_rank = median(rank_pearson),
  q25_rank = quantile(rank_pearson, 0.25),
  q75_rank = quantile(rank_pearson, 0.75),
  worst_rank = max(rank_pearson),
  best_rank = min(rank_pearson),
  top10_frequency = mean(rank_pearson <= 10),
  top25_frequency = mean(rank_pearson <= 25),
  median_pearson = median(pearson)
), by = ligand][order(median_rank)]
fwrite(loo_summary, file.path(out_dir, "nichenet_leave_one_donor_out_summary.csv"))

# Figure 1: primary ligand activity ranking.
plot_primary <- head(primary, 25)
plot_primary$ligand <- factor(plot_primary$ligand, levels = rev(plot_primary$ligand))
plot_primary$highlight <- ifelse(plot_primary$ligand %in% c("SPP1", "S100A8", "S100A9", "TGFB1", "OSM", "LIF", "AREG"),
                                 "Prespecified", "Other")
p1 <- ggplot(plot_primary, aes(x = pearson, y = ligand, fill = highlight)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c(Prespecified = "#D62728", Other = "#4C78A8")) +
  labs(title = "NicheNet ligand activity for the IPF secretory epithelial program",
       subtitle = "Full-transcriptome paired Club/Goblet vs ATII program; 15 IPF donors",
       x = "Pearson ligand activity", y = NULL, fill = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(out_dir, "Figure_1_NicheNet_ligand_activity.png"), p1, width = 8.5, height = 7.5, dpi = 400)
ggsave(file.path(out_dir, "Figure_1_NicheNet_ligand_activity.pdf"), p1, width = 8.5, height = 7.5)

# Figure 2: regulatory potential for axis genes (column-scaled for visibility).
heat_focus <- unique(c("SPP1", head(primary$ligand, 7)))
heat_focus <- intersect(heat_focus, focus_ligands)
heat_dt <- axis_links[ligand %in% heat_focus]
heat_dt[, scaled_weight := {
  rng <- range(regulatory_potential, na.rm = TRUE)
  if (diff(rng) == 0) rep(0, .N) else (regulatory_potential - rng[1]) / diff(rng)
}, by = ligand]
p2 <- ggplot(heat_dt, aes(x = ligand, y = target, fill = scaled_weight)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = ifelse(global_target_rank <= 250, paste0("#", global_target_rank), "")), size = 2.6) +
  scale_fill_gradient(low = "white", high = "#B2182B", na.value = "grey95") +
  labs(title = "NicheNet prior support for epithelial-axis targets",
       subtitle = "Numbers mark targets ranked within the ligand's global top 250",
       x = "Candidate ligand", y = NULL, fill = "Within-ligand\nscaled weight") +
  theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 40, hjust = 1))
ggsave(file.path(out_dir, "Figure_2_NicheNet_axis_target_heatmap.png"), p2, width = 9.5, height = 6.8, dpi = 400)
ggsave(file.path(out_dir, "Figure_2_NicheNet_axis_target_heatmap.pdf"), p2, width = 9.5, height = 6.8)

# Figure 3: SPP1 and other prespecified ligands across geneset definitions.
p3_dt <- sensitivity[geneset %in% c("Secretory_up_top250", "Secretory_up_top500_primary",
                                     "Secretory_up_all_FDR05_logFC1", "ATII_down_top500")]
p3 <- ggplot(p3_dt, aes(x = geneset, y = rank_pearson, group = ligand, color = ligand)) +
  geom_line(linewidth = 0.7, alpha = 0.8) + geom_point(size = 2.4) +
  scale_y_reverse() +
  labs(title = "Sensitivity of prespecified ligand ranks",
       subtitle = "Lower rank is better; NicheNet regulatory potential is non-directional",
       x = NULL, y = "Pearson activity rank", color = "Ligand") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "right")
ggsave(file.path(out_dir, "Figure_3_NicheNet_rank_sensitivity.png"), p3, width = 10.5, height = 6.5, dpi = 400)
ggsave(file.path(out_dir, "Figure_3_NicheNet_rank_sensitivity.pdf"), p3, width = 10.5, height = 6.5)

# Figure 4: donor robustness for the most stable and prespecified ligands.
stable_ligands <- head(loo_summary$ligand, 10)
robust_ligands <- unique(c(stable_ligands, intersect(prespecified, loo_summary$ligand)))
p4_dt <- loo_activity[ligand %in% robust_ligands]
rank_order <- loo_summary[ligand %in% robust_ligands][order(median_rank), ligand]
p4_dt$ligand <- factor(p4_dt$ligand, levels = rev(rank_order))
p4 <- ggplot(p4_dt, aes(x = rank_pearson, y = ligand)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, fill = "#D9EAF7") +
  geom_jitter(width = 0, height = 0.14, size = 1.6, alpha = 0.65, color = "#2C5D8A") +
  geom_vline(xintercept = 25.5, linetype = "dashed", color = "grey45") +
  labs(title = "Leave-one-donor-out robustness of NicheNet ligand ranks",
       subtitle = "Each point drops one of 15 paired IPF donors",
       x = "Pearson activity rank (lower is better)", y = NULL) +
  theme_bw(base_size = 11)
ggsave(file.path(out_dir, "Figure_4_NicheNet_leave_one_donor_out.png"), p4, width = 9, height = 7.5, dpi = 400)
ggsave(file.path(out_dir, "Figure_4_NicheNet_leave_one_donor_out.pdf"), p4, width = 9, height = 7.5)

# Figure 5: integrate sender expression, receiver receptor support and ligand activity.
p5_dt <- copy(focused_source)
p5_dt$ligand <- factor(p5_dt$ligand, levels = p5_dt[order(pearson), ligand])
p5 <- ggplot(p5_dt, aes(x = ligand, y = pearson,
                        size = ligand_detect_fraction_myeloid,
                        fill = myeloid_vs_secretory_detection_log2ratio)) +
  geom_point(shape = 21, color = "black", stroke = 0.5) +
  geom_text(aes(label = paste0("rank ", rank_pearson)), vjust = -1.35, size = 3.2) +
  scale_size_continuous(range = c(5, 12), limits = c(0, 1), name = "Myeloid\ndetection") +
  scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0,
                       name = "log2 myeloid /\nsecretory detection") +
  labs(title = "Integrated evidence for prespecified upstream ligands",
       subtitle = "Activity predicts the secretory program; size and color reflect sender localization",
       x = NULL, y = "NicheNet Pearson ligand activity") +
  theme_bw(base_size = 11) + theme(legend.position = "right") +
  expand_limits(y = max(p5_dt$pearson) * 1.18)
ggsave(file.path(out_dir, "Figure_5_NicheNet_focused_ligand_evidence.png"), p5, width = 9.5, height = 6.5, dpi = 400)
ggsave(file.path(out_dir, "Figure_5_NicheNet_focused_ligand_evidence.pdf"), p5, width = 9.5, height = 6.5)

spp1_primary <- primary[ligand == "SPP1"]
writeLines(c(
  "NicheNet GSE136831 sender-focused analysis",
  "NicheNet human prior source: https://zenodo.org/records/7074291",
  sprintf("Ligand-target matrix MD5: %s", unname(tools::md5sum(file.path(data_dir, "ligand_target_matrix_nsga2r_final.rds")))),
  sprintf("Ligand-receptor network MD5: %s", unname(tools::md5sum(file.path(data_dir, "lr_network_human_21122021.rds")))),
  sprintf("Paired IPF donors: %d", length(paired_donors)),
  sprintf("Receiver background genes in NicheNet model: %d", length(background)),
  sprintf("Primary secretory-up genes: %d", length(geneset_up_500)),
  sprintf("Expression-supported potential ligands: %d", length(potential_ligands)),
  sprintf("SPP1 primary Pearson rank: %s/%d", ifelse(nrow(spp1_primary), spp1_primary$rank_pearson, "not eligible"), length(potential_ligands)),
  sprintf("SPP1 primary Pearson activity: %s", ifelse(nrow(spp1_primary), signif(spp1_primary$pearson, 4), "NA")),
  sprintf("SPP1 top-250 target overlap: %s; enrichment FDR: %s",
          ifelse(nrow(spp1_primary), spp1_primary$top250_overlap, "NA"),
          ifelse(nrow(spp1_primary), signif(spp1_primary$fisher_FDR, 4), "NA")),
  sprintf("SPP1 leave-one-donor-out median rank: %s",
          ifelse("SPP1" %in% loo_summary$ligand, loo_summary[ligand == "SPP1", median_rank], "NA")),
  sprintf("TGFB1 leave-one-donor-out median rank: %s; myeloid detection: %.3f; secretory detection: %.3f",
          loo_summary[ligand == "TGFB1", median_rank], ligand_detect_fraction["TGFB1"],
          ligand_detect_fraction_secretory["TGFB1"]),
  sprintf("S100A8 leave-one-donor-out median rank: %s; myeloid detection: %.3f; secretory detection: %.3f",
          loo_summary[ligand == "S100A8", median_rank], ligand_detect_fraction["S100A8"],
          ligand_detect_fraction_secretory["S100A8"]),
  "Interpretation: NicheNet scores whether a ligand's prior target profile predicts the observed receiver gene set.",
  "The model is hypothesis-generating and its regulatory-potential weights are not directional causal effects."
), file.path(out_dir, "README.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_sessionInfo.txt"))

cat("Receiver DE: up=", length(up_sig), " down=", length(down_sig), " background=", length(background), "\n", sep = "")
cat("Potential ligands:", length(potential_ligands), "\n")
cat("\nPrimary top 25:\n")
print(primary[1:min(25, .N), .(ligand, pearson, aupr_corrected, top250_overlap, fisher_FDR,
                                rank_pearson, ligand_detect_fraction_myeloid, n_expressed_receptors)])
cat("\nPrespecified sensitivity:\n")
print(sensitivity[order(geneset, rank_pearson)])
cat("\nSPP1 axis links:\n")
print(axis_links[ligand == "SPP1"][order(global_target_rank)])
cat("\nLeave-one-donor-out summary for prespecified ligands:\n")
print(loo_summary[ligand %in% prespecified][order(median_rank)])

