# Background deletion QC (plan §6.6, D12): only meaningful on counts made with count_all_bases.
# Site calling always filters ref == "T", so this never changes the calls.

summarise_delrate <- function(dt, by) {
  dt[, .(n_sites = .N, median = stats::median(delrate), q1 = stats::quantile(delrate, 0.25, names = FALSE),
         q3 = stats::quantile(delrate, 0.75, names = FALSE), iqr = stats::IQR(delrate), mean = mean(delrate)),
     keyby = by]
}

#' Background deletion-rate summary per sample and reference base, and input vs BS.
#'
#' @param merged counts_merged.tsv path or data.table (needs sample_id, ref, delrate; treat optional).
#' @param outdir Directory for background_qc.tsv, background_qc_by_treat.tsv and
#'   background_qc.pdf/.png (boxplot faceted by reference base).
#' @return list of the two summary tables (invisibly).
background_qc <- function(merged, outdir = "background") {
  dt <- if (is.character(merged)) read_text_table(merged) else copy(as.data.table(merged))
  for (c in c("sample_id", "ref", "delrate")) if (!c %in% names(dt)) stop("merged counts lack column ", c, call. = FALSE)
  dt[, delrate := as.numeric(delrate)]
  if (all(dt$ref %in% "T")) warning("only T sites present; run counting with count_all_bases for background QC", call. = FALSE)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  per_sample <- summarise_delrate(dt, c("sample_id", "ref"))
  fwrite(per_sample, file.path(outdir, "background_qc.tsv"), sep = "\t", quote = FALSE, na = "NA")
  by_treat <- NULL
  if ("treat" %in% names(dt)) {
    by_treat <- summarise_delrate(dt, c("ref", "treat"))
  } else {
    warning("no treat column; input vs BS summary skipped", call. = FALSE)
    by_treat <- dt[0, .(ref, treat = character())]
  }
  fwrite(by_treat, file.path(outdir, "background_qc_by_treat.tsv"), sep = "\t", quote = FALSE, na = "NA")

  plot_dt <- copy(dt)
  if (!"treat" %in% names(plot_dt)) plot_dt[, treat := "all"]
  p <- ggplot2::ggplot(plot_dt, ggplot2::aes(x = sample_id, y = delrate, fill = treat)) +
    ggplot2::geom_boxplot(outlier.size = 0.5) +
    ggplot2::facet_wrap(~ref, scales = "free_y") +
    ggplot2::labs(x = "Sample", y = "Deletion rate", fill = "Treatment",
                  title = "Background deletion rate by reference base") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
  w <- max(7, 0.25 * length(unique(dt$sample_id)) + 3)
  ggplot2::ggsave(file.path(outdir, "background_qc.pdf"), p, width = w, height = 6)
  ggplot2::ggsave(file.path(outdir, "background_qc.png"), p, width = w, height = 6, dpi = 150)
  invisible(list(per_sample = per_sample, by_treat = by_treat))
}
