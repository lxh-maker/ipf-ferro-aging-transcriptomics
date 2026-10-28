suppressPackageStartupMessages({
  library(Matrix)
  library(scTenifoldNet)
  library(scTenifoldKnk)
  library(ggplot2)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
root <- normalizePath(file.path(dirname(script_file), "..", ".."), mustWork = FALSE)
ko_root <- file.path(root, "03_virtual_knockout", "results")
out_root <- file.path(root, "04_virtual_overexpression", "results")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

jobs <- list(
  list(compartment = "IPF_Macrophages", target = "S100A8"),
  list(compartment = "IPF_Epithelial", target = "LCN2"),
  list(compartment = "IPF_Epithelial", target = "ACSL4"),
  list(compartment = "IPF_Fibroblast", target = "LIFR")
)
doses <- c(1.5, 2, 3)
n_cores <- max(1L, min(8L, parallel::detectCores() - 2L))
summary_rows <- list()

for (job in jobs) {
  compartment <- job$compartment
  target <- job$target
  ko_dir <- file.path(ko_root, paste0(compartment, "__", target, "_KO"))
  ko <- readRDS(file.path(ko_dir, "scTenifoldKnk_result.rds"))
  WT <- as.matrix(ko$tensorNetworks$WT)
  stopifnot(target %in% rownames(WT))

  # scTenifoldKnk infers KO directions from local expression structure.  The
  # sign below is the first-order inverse used only to annotate an exploratory
  # OE sensitivity analysis; the manifold distances are calculated directly
  # from amplified outgoing edges and do not depend on this sign heuristic.
  ko_tab <- ko$diffRegulation
  inverse_direction <- setNames(-ko_tab$directionScore, ko_tab$gene)

  job_dir <- file.path(out_root, paste0(compartment, "__", target, "_OE"))
  dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
  all_doses <- list()

  for (dose in doses) {
    message("Running ", compartment, " / ", target, " OE ", dose, "x")
    OE <- WT
    OE[target, ] <- WT[target, ] * dose
    set.seed(20261027)
    MA <- scTenifoldNet::manifoldAlignment(WT, OE, d = 2, nCores = n_cores)
    dose_direction <- inverse_direction * (dose - 1)
    DR <- scTenifoldKnk::dRegulation(
      MA,
      gKO = target,
      empiricalNull = FALSE,
      direction = dose_direction
    )
    DR <- as.data.frame(DR)
    DR <- DR[DR$gene != target, , drop = FALSE]
    DR$dose <- dose
    DR <- DR[order(DR$p.adj, -DR$distance), , drop = FALSE]
    all_doses[[as.character(dose)]] <- DR
    write.csv(
      DR,
      file.path(job_dir, paste0("OE_", dose, "x_differential_regulation.csv")),
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      compartment = compartment,
      target = target,
      dose = dose,
      significant_responses_FDR05 = sum(DR$p.adj < 0.05, na.rm = TRUE),
      top_response = if (nrow(DR)) DR$gene[1] else NA_character_,
      top_FDR = if (nrow(DR)) DR$p.adj[1] else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  combined <- do.call(rbind, all_doses)
  rownames(combined) <- NULL
  write.csv(combined, file.path(job_dir, "all_doses.csv"), row.names = FALSE, fileEncoding = "UTF-8")

  stable <- Reduce(intersect, lapply(all_doses, function(z) z$gene[z$p.adj < 0.05]))
  stable_tab <- combined[combined$gene %in% stable, , drop = FALSE]
  write.csv(stable_tab, file.path(job_dir, "stable_responses_FDR05_all_doses.csv"), row.names = FALSE, fileEncoding = "UTF-8")

  top_genes <- unique(unlist(lapply(all_doses, function(z) head(z$gene, 12))))
  plot_tab <- combined[combined$gene %in% top_genes, , drop = FALSE]
  plot_tab$score <- -log10(pmax(plot_tab$p.adj, .Machine$double.xmin))
  plot_tab$gene <- factor(plot_tab$gene, levels = rev(unique(top_genes)))
  p <- ggplot(plot_tab, aes(x = factor(dose), y = gene, size = score, color = directionScore)) +
    geom_point(alpha = 0.85) +
    scale_color_gradient2(low = "#2563EB", mid = "grey90", high = "#DC2626", midpoint = 0) +
    scale_size_continuous(range = c(1, 8)) +
    theme_bw(base_size = 11) +
    labs(
      title = paste0(target, " virtual overexpression sensitivity in ", compartment),
      subtitle = "Outgoing-edge amplification; direction is an inverse-KO first-order estimate",
      x = "Network amplification factor",
      y = NULL,
      size = "-log10(FDR)",
      color = "Estimated direction"
    )
  ggsave(file.path(job_dir, "dose_response.pdf"), p, width = 8, height = 7)
  ggsave(file.path(job_dir, "dose_response.png"), p, width = 8, height = 7, dpi = 300)
  writeLines(capture.output(sessionInfo()), file.path(job_dir, "sessionInfo.txt"), useBytes = TRUE)
}

summary_table <- do.call(rbind, summary_rows)
write.csv(summary_table, file.path(out_root, "virtual_overexpression_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")
print(summary_table)

