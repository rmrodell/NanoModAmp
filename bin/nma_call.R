#!/usr/bin/env Rscript
# nma_call.R: run one site-calling analysis from the analyses YAML (plan §6.3, §6.4, §5.4, §5.5).
# Writes <outdir>/<name>/{<name>_log.txt, <name>_resolved_config.yaml, data_summary/, data_raw/, plots/}.
suppressPackageStartupMessages({
  library(argparse)
  library(nanomodamp)
})

p <- ArgumentParser(description = "Run one NanoModAmp site-calling analysis.")
p$add_argument("--counts", required = TRUE, help = "counts_merged.tsv (§5.3)")
p$add_argument("--analyses", required = TRUE, help = "analyses YAML (§5.4)")
p$add_argument("--name", required = TRUE, help = "name of the analysis to run")
p$add_argument("--outdir", default = ".", help = "directory that receives <name>/ [default: .]")
p$add_argument("--cpus", type = "integer", default = 1L, help = "parallel workers (task.cpus) [default: 1]")
p$add_argument("--p_adjust", default = "BH", choices = c("BH", "legacy"),
               help = "default for analyses that do not set p_adjust (R-30) [default: BH]")
p$add_argument("--plot_all_sites", action = "store_true", default = FALSE,
               help = "default for analyses that do not set plot_all_sites")
args <- p$parse_args()

cfg <- read_analyses(args$analyses, plot_all_sites = args$plot_all_sites, p_adjust = args$p_adjust)
if (!args$name %in% names(cfg$analyses)) {
  stop("analysis '", args$name, "' is not in ", args$analyses, call. = FALSE)
}
counts <- read_counts(args$counts)
run_analysis(counts, cfg$analyses[[args$name]], args$outdir, workers = max(1L, args$cpus))
cat("nma_call.R: finished", args$name, "\n")
