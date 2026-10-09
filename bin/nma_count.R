#!/usr/bin/env Rscript
# Per-sample site counts (plan §6.2): writes <prefix>.counts.tsv and <prefix>.failed_regions.tsv.
suppressPackageStartupMessages({ library(optparse); library(nanomodamp) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--bam", help = "Indexed BAM"),
  make_option("--fasta", help = "Indexed reference FASTA"),
  make_option("--bed", help = "BED (>= 6 columns)"),
  make_option("--prefix", help = "Output prefix (sample_id)"),
  make_option("--bed-coordinates", dest = "bed_coordinates", default = "bed0", help = "bed0 | one_based_start [%default]"),
  make_option("--min-coverage", dest = "min_coverage", type = "integer", default = 20L, help = "keep totalReads > this [%default]"),
  make_option("--max-depth", dest = "max_depth", type = "integer", default = 200000L, help = "[%default]"),
  make_option("--min-mapq", dest = "min_mapq", type = "integer", default = 1L, help = "[%default]"),
  make_option("--min-base-quality", dest = "min_base_quality", type = "integer", default = 1L, help = "[%default]"),
  make_option("--all-bases", dest = "all_bases", action = "store_true", default = FALSE, help = "report all reference bases (D12)"),
  make_option("--allow-region-failures", dest = "allow_region_failures", action = "store_true", default = FALSE),
  make_option("--threads", type = "integer", default = 1L, help = "[%default]")
)))
for (o in c("bam", "fasta", "bed", "prefix")) if (is.null(opt[[o]])) stop("--", o, " is required", call. = FALSE)

failed_file <- paste0(opt$prefix, ".failed_regions.tsv")
res <- tryCatch(
  count_sites(opt$bam, opt$fasta, opt$bed, bed_coordinates = opt$bed_coordinates,
              min_coverage = opt$min_coverage, max_depth = opt$max_depth, min_mapq = opt$min_mapq,
              min_base_quality = opt$min_base_quality, count_all_bases = opt$all_bases,
              threads = opt$threads, allow_region_failures = opt$allow_region_failures),
  nma_region_failure = function(e) {
    data.table::fwrite(e$failed, failed_file, sep = "\t", quote = FALSE)
    message("ERROR: ", conditionMessage(e), " (see ", failed_file, "; --allow-region-failures to continue)")
    quit(status = 1)
  })
write_counts(res$counts, paste0(opt$prefix, ".counts.tsv"))
data.table::fwrite(res$failed, failed_file, sep = "\t", quote = FALSE)
if (nrow(res$failed)) message("WARNING: ", nrow(res$failed), " region(s) failed; listed in ", failed_file)
if (res$duplicate_sites > 0) message("NOTE: ", res$duplicate_sites, " site row(s) repeated because BED regions overlap (kept as in legacy)")
if (nrow(res$saturated)) message("WARNING: pileup max_depth (", opt$max_depth, ") reached at ", nrow(res$saturated),
                                 " position(s), e.g. ", res$saturated$chr[1], ":", res$saturated$pos[1], " (R-16)")
message(opt$prefix, ": ", nrow(res$counts), " site rows")
