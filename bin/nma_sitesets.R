#!/usr/bin/env Rscript
# nma_sitesets.R: build the site sets of an analyses YAML from finished calling/<name>/ results
# (plan §6.5). Writes <outdir>/<set>/{data_summary/, data_raw/, plots/}.
suppressPackageStartupMessages({
  library(argparse)
  library(nanomodamp)
})

p <- ArgumentParser(description = "Build NanoModAmp site sets (union / intersection / difference).")
p$add_argument("--analyses", required = TRUE, help = "analyses YAML (§5.4)")
p$add_argument("--calling_dir", required = TRUE, help = "directory holding the calling/<name>/ results")
p$add_argument("--counts", default = NULL, help = "counts_merged.tsv, for the raw-data tables (optional)")
p$add_argument("--outdir", default = "site_sets", help = "output directory [default: site_sets]")
args <- p$parse_args()

cfg <- read_analyses(args$analyses)
counts <- if (is.null(args$counts)) NULL else read_counts(args$counts)
run_site_sets(cfg, args$calling_dir, args$outdir, counts)
cat("nma_sitesets.R: built", length(cfg$site_sets), "site set(s)\n")
