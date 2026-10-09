#!/usr/bin/env Rscript
# Background deletion QC (plan §6.6; count_all_bases only).
#   nma_bgqc.R counts_merged.tsv [outdir]
suppressPackageStartupMessages({ library(nanomodamp) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: nma_bgqc.R counts_merged.tsv [outdir]", call. = FALSE)
background_qc(args[1], if (length(args) > 1) args[2] else "background")
