#!/usr/bin/env Rscript
# Merge count tables by column name (plan §5.3, R-06).
#   per-sample mode:  nma_merge.R --samplesheet samplesheet.csv --counts a.counts.tsv b.counts.tsv ...
#   run-merge mode (D32, --input_counts):  nma_merge.R --tables t1.tsv t2.tsv ... [--sources s1 s2 ...]
suppressPackageStartupMessages({ library(nanomodamp) })

args <- commandArgs(trailingOnly = TRUE)
take <- function(flag) {
  i <- match(flag, args)
  if (is.na(i)) return(NULL)
  j <- i + 1
  out <- character()
  while (j <= length(args) && !startsWith(args[j], "--")) { out <- c(out, args[j]); j <- j + 1 }
  out
}
out_file <- c(take("--out"), "counts_merged.tsv")[1]
if (!is.null(take("--tables"))) {
  tables <- take("--tables")
  sources <- take("--sources"); if (is.null(sources)) sources <- tables
  res <- merge_count_tables(tables, sources)
  write_counts(res$merged, out_file)
  data.table::fwrite(res$sources, c(take("--sources-out"), "merge_sources.tsv")[1], sep = "\t", quote = FALSE)
  message("merged ", length(tables), " tables: ", nrow(res$sources), " samples, ", nrow(res$merged), " rows")
} else {
  ss <- take("--samplesheet"); counts <- take("--counts")
  if (is.null(ss) || is.null(counts)) stop("need --samplesheet and --counts, or --tables", call. = FALSE)
  merged <- merge_counts(counts, ss)
  write_counts(merged, out_file)
  message("merged ", length(counts), " samples: ", nrow(merged), " rows")
}
