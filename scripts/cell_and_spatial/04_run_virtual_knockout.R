suppressPackageStartupMessages({
  library(Matrix)
  library(scTenifoldKnk)
  library(ggplot2)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
root <- normalizePath(file.path(dirname(script_file), "..", ".."), mustWork = FALSE)
input_root <- file.path(root, "03_virtual_knockout", "inputs")
output_root <- file.path(root, "03_virtual_knockout", "results")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

jobs <- list(
  list(compartment = "IPF_Macrophages", target = "S100A8"),
  list(compartment = "IPF_Epithelial", target = "LCN2"),
  list(compartment = "IPF_Epithelial", target = "ACSL4"),
  list(compartment = "IPF_Fibroblast", target = "LIFR")
)

read_input <- function(compartment) {
  path <- file.path(input_root, compartment)
  counts <- Matrix::readMM(file.path(path, "counts.mtx"))
  genes <- readLines(file.path(path, "genes.tsv"), encoding = "UTF-8", warn = FALSE)
  cells <- readLines(file.path(path, "cells.tsv"), encoding = "UTF-8", warn = FALSE)
  rownames(counts) <- genes
  colnames(counts) <- cells
  keep <- Matrix::rowSums(counts) > 0
  counts <- counts[keep, , drop = FALSE]
  as.matrix(counts)
}

save_knockout_plots <- function(tab, target, compartment, out_dir) {
  tab$minus_log10_fdr <- -log10(pmax(tab$p.adj, .Machine$double.xmin))
  tab$significant <- tab$p.adj < 0.05
  xcol <- if ("directionScore" %in% colnames(tab)) "directionScore" else if ("Z" %in% colnames(tab)) "Z" else "distance"
  p <- ggplot(tab, aes(x = .data[[xcol]], y = minus_log10_fdr, color = significant)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_manual(values = c(`FALSE` = "grey65", `TRUE` = "#D62728")) +
    theme_bw(base_size = 12) +
    labs(
      title = paste0(target, " virtual knockout in ", compartment),
      x = xcol,
      y = "-log10(FDR)",
      color = "FDR < 0.05"
    )
  ggsave(file.path(out_dir, "knockout_scatter.pdf"), p, width = 7, height = 5)
  ggsave(file.path(out_dir, "knockout_scatter.png"), p, width = 7, height = 5, dpi = 300)

  # Display the strongest statistically supported responses. Ranking only by
  # direction magnitude can promote non-significant genes into this panel.
  direction_magnitude <- if ("directionScore" %in% colnames(tab)) abs(tab$directionScore) else if ("Z" %in% colnames(tab)) abs(tab$Z) else tab$distance
  top <- tab[order(tab$p.adj, -direction_magnitude), , drop = FALSE]
  top <- head(top, 20)
  top$gene <- factor(top$gene, levels = rev(top$gene))
  top$value <- if ("directionScore" %in% colnames(top)) top$directionScore else if ("Z" %in% colnames(top)) top$Z else top$distance
  p2 <- ggplot(top, aes(x = gene, y = value, fill = value > 0)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c(`FALSE` = "#3B82F6", `TRUE` = "#DC2626"), guide = "none") +
    theme_bw(base_size = 11) +
    labs(title = paste0("Top predicted responses to ", target, " knockout"), x = NULL, y = xcol)
  ggsave(file.path(out_dir, "top20_responses.pdf"), p2, width = 7, height = 6)
  ggsave(file.path(out_dir, "top20_responses.png"), p2, width = 7, height = 6, dpi = 300)
}

summary_rows <- list()
n_cores <- max(1L, min(8L, parallel::detectCores() - 2L))

for (job in jobs) {
  compartment <- job$compartment
  target <- job$target
  out_dir <- file.path(output_root, paste0(compartment, "__", target, "_KO"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  result_file <- file.path(out_dir, "scTenifoldKnk_result.rds")
  table_file <- file.path(out_dir, "differential_regulation.csv")

  message("Running ", compartment, " / ", target)
  counts <- read_input(compartment)
  if (!target %in% rownames(counts)) {
    warning(target, " is absent from ", compartment)
    next
  }

  start_time <- Sys.time()
  result <- tryCatch(
    scTenifoldKnk(
      countMatrix = counts,
      gKO = target,
      transcriptomeWide = FALSE,
      qc = FALSE,
      nc_nNet = 10,
      nc_nCells = min(500L, ncol(counts)),
      nc_nComp = 3,
      nc_q = 0.90,
      td_K = 3,
      ma_nDim = 2,
      ma_method = "manifold",
      dr_empiricalNull = FALSE,
      dr_direction = TRUE,
      nCores = n_cores,
      seed = 20261027
    ),
    error = function(e) {
      writeLines(conditionMessage(e), file.path(out_dir, "ERROR.txt"), useBytes = TRUE)
      NULL
    }
  )
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  if (is.null(result)) next

  saveRDS(result, result_file)
  tab <- as.data.frame(result$diffRegulation)
  if (!"gene" %in% colnames(tab)) tab$gene <- rownames(tab)
  tab <- tab[tab$gene != target, , drop = FALSE]
  tab <- tab[order(tab$p.adj, -tab$distance), , drop = FALSE]
  write.csv(tab, table_file, row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(tab[tab$p.adj < 0.05, , drop = FALSE], file.path(out_dir, "significant_responses_FDR05.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  save_knockout_plots(tab, target, compartment, out_dir)

  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    compartment = compartment,
    target = target,
    cells = ncol(counts),
    genes = nrow(counts),
    target_detected_cells = sum(counts[target, ] > 0),
    significant_responses_FDR05 = sum(tab$p.adj < 0.05, na.rm = TRUE),
    top_response = if (nrow(tab)) tab$gene[1] else NA_character_,
    top_FDR = if (nrow(tab)) tab$p.adj[1] else NA_real_,
    runtime_minutes = elapsed,
    stringsAsFactors = FALSE
  )
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"), useBytes = TRUE)
}

if (length(summary_rows)) {
  summary_table <- do.call(rbind, summary_rows)
  write.csv(summary_table, file.path(output_root, "virtual_knockout_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  print(summary_table)
}
