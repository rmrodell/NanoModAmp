# Per-site counting (plan §6.2), ported from PUS7regulation2026 Figure2/BIDdetect/bam_counts_fast.R.
#
# Kept exactly as written: pileup parameters, strand filter and complementing, wide reshape,
# ref/kmer lookup, T-only filter, totalReads = A+C+G+T+Deletion+Insertion (D9), totalReads >
# min_coverage (D10), delrate = Deletion/totalReads, row order (BED order, then position).
# Approved changes: BED coordinate convention is a parameter (R-15/D27, default bed0); region
# errors are collected and fail the task unless allowed (R-12); max_depth saturation warns (R-16);
# output always has the contract header (§5.3), `ref` column name, kmer written as NA.

#' Count columns of a per-sample counts table, in contract order (§5.3).
COUNT_COLUMNS <- c("chr", "pos", "gene", "totalReads", "A.count", "C.count", "G.count",
                   "T.count", "Deletion.count", "Insertion.count", "ref", "kmer", "strand",
                   "delrate")

BASE_COLUMNS <- c("A", "C", "G", "T", "Deletion", "Insertion")

#' Read a BED file into regions with 1-based, inclusive starts.
#'
#' @param bed Path to a BED file with at least 6 columns (chr, start, end, name, score, strand).
#' @param bed_coordinates `"bed0"` (standard, default): start is 0-based, so 1 is added and
#'   rows with start = end are an error. `"one_based_start"`: legacy reading, start and end are
#'   used as written (1-based, inclusive), as in the paper's single-site BEDs (R-15).
#' @return data.table with chr, start (1-based), end, gene, strand, bed_start, bed_end.
read_bed_regions <- function(bed, bed_coordinates = c("bed0", "one_based_start")) {
  bed_coordinates <- match.arg(bed_coordinates)
  b <- fread(bed, header = FALSE, sep = "\t", quote = "", colClasses = list(character = c(1, 4, 6)))
  if (ncol(b) < 6) stop("BED file has fewer than 6 columns: ", bed, call. = FALSE)
  b <- b[, .(chr = V1, bed_start = as.integer(V2), bed_end = as.integer(V3), gene = V4, strand = V6)]
  if (any(b$bed_end < b$bed_start)) stop("BED rows with end < start in ", bed, call. = FALSE)
  if (bed_coordinates == "bed0") {
    zero <- b[bed_start == bed_end]
    if (nrow(zero) > 0) {
      stop(sprintf(paste0("%d BED row(s) have start = end (e.g. %s:%d), which is empty in standard BED ",
                          "(bed_coordinates = bed0). For a legacy single-site BED use ",
                          "bed_coordinates = one_based_start (R-15)."),
                   nrow(zero), zero$chr[1], zero$bed_start[1]), call. = FALSE)
    }
    b[, start := bed_start + 1L]
  } else {
    b[, start := bed_start]
  }
  b[, end := bed_end]
  b[, .(chr, start, end, gene, strand, bed_start, bed_end)]
}

# One region, as process_region() in bam_counts_fast.R.
count_region <- function(region, fasta_handle, bam_handle, all_bases, pileup_params, max_depth) {
  region_gr <- GRanges(seqnames = region$chr,
                       ranges = IRanges(start = region$start, end = region$end),
                       strand = region$strand)
  PU <- as.data.table(pileup(bam_handle, scanBamParam = ScanBamParam(which = region_gr),
                             pileupParam = pileup_params))
  if (nrow(PU) == 0) return(NULL)
  if (region$strand != "*") PU <- PU[strand == region$strand]
  if (nrow(PU) == 0) return(NULL)

  PU[, nucleotide_new := as.character(nucleotide)]
  PU[strand == "-", nucleotide_new := chartr("ATCG", "TAGC", nucleotide_new)]
  # R-16: reads stacked at a position, excluding insertions (pileup's max_depth caps these)
  depth <- PU[nucleotide_new != "+", .(n = sum(count)), by = pos]
  saturated <- depth[n >= max_depth, pos]

  counts_dt <- dcast(PU, seqnames + pos ~ nucleotide_new, value.var = "count",
                     fun.aggregate = sum, fill = 0L)
  setnames(counts_dt, "seqnames", "chr")
  if ("-" %in% names(counts_dt)) setnames(counts_dt, "-", "Deletion")
  if ("+" %in% names(counts_dt)) setnames(counts_dt, "+", "Insertion")

  pos_dt <- data.table(pos = unique(counts_dt$pos))
  ref_seqs <- getSeq(fasta_handle, GRanges(seqnames = region$chr,
                                           ranges = IRanges(start = pos_dt$pos, end = pos_dt$pos)))
  if (region$strand == "-") ref_seqs <- complement(ref_seqs)
  pos_dt[, ref := as.character(ref_seqs)]

  pos_dt[, kmer := NA_character_]
  chrom_len <- seqlengths(fasta_handle)[region$chr]
  ok <- (pos_dt$pos >= 3) & (pos_dt$pos <= (chrom_len - 2))
  if (any(ok)) {
    p <- pos_dt$pos[ok]
    kmer_seqs <- getSeq(fasta_handle, GRanges(seqnames = region$chr,
                                              ranges = IRanges(start = p - 2, end = p + 2)))
    if (region$strand == "-") kmer_seqs <- reverseComplement(kmer_seqs)
    pos_dt[ok, kmer := as.character(kmer_seqs)]
  }

  res <- merge(counts_dt, pos_dt, by = "pos", all.x = TRUE)
  if (!all_bases) res <- res[ref == "T"]
  if (nrow(res) == 0) return(list(rows = NULL, saturated = saturated))
  res[, gene := region$gene]
  res[, strand := region$strand]
  res[, chr := region$chr]
  list(rows = res, saturated = saturated)
}

#' Count bases, deletions and insertions per site (plan §6.2).
#'
#' @param bam Indexed BAM file.
#' @param fasta Indexed FASTA (`.fai` next to it).
#' @param bed BED file path, or regions from [read_bed_regions()].
#' @param bed_coordinates See [read_bed_regions()].
#' @param min_coverage Keep sites with `totalReads > min_coverage` (D10).
#' @param max_depth,min_mapq,min_base_quality Pileup parameters (paper: 200000, 1, 1).
#' @param count_all_bases Report all reference bases, not only T (D12).
#' @param threads Parallel workers over regions; results do not depend on it.
#' @param allow_region_failures If FALSE, any failed region is an error (R-12).
#' @return list: `counts` (data.table, §5.3 columns), `failed` (chr, start, end, gene, error, in
#'   BED coordinates as given), `duplicate_sites` (number of (chr, pos) rows from overlapping
#'   regions), `saturated` (data.table of chr, pos at max_depth).
count_sites <- function(bam, fasta, bed, bed_coordinates = "bed0", min_coverage = 20,
                        max_depth = 200000L, min_mapq = 1L, min_base_quality = 1L,
                        count_all_bases = FALSE, threads = 1L, allow_region_failures = FALSE) {
  regions <- if (is.character(bed)) read_bed_regions(bed, bed_coordinates) else as.data.table(bed)
  fasta_handle <- FaFile(fasta)
  bam_handle <- BamFile(bam)
  pp <- PileupParam(max_depth = as.integer(max_depth), min_mapq = as.integer(min_mapq),
                    distinguish_nucleotides = TRUE, ignore_query_Ns = TRUE,
                    min_base_quality = as.integer(min_base_quality),
                    include_insertions = TRUE, include_deletions = TRUE, distinguish_strands = TRUE)

  one <- function(i) {
    r <- regions[i]
    tryCatch(c(count_region(r, fasta_handle, bam_handle, count_all_bases, pp, max_depth), list(error = NA_character_)),
             error = function(e) list(rows = NULL, saturated = integer(0), error = conditionMessage(e)))
  }
  idx <- seq_len(nrow(regions))
  if (threads > 1 && length(idx) > 1) {
    old <- future::plan(future::multisession, workers = as.integer(threads))
    on.exit(future::plan(old), add = TRUE)
    res <- future.apply::future_lapply(idx, one, future.seed = TRUE)
  } else {
    res <- lapply(idx, one)
  }

  err <- vapply(res, function(x) if (is.null(x$error)) NA_character_ else x$error, character(1))
  failed <- regions[!is.na(err), .(chr, start = bed_start, end = bed_end, gene)]
  failed[, error := err[!is.na(err)]]
  if (nrow(failed) > 0 && !allow_region_failures) {
    msg <- sprintf("%d region(s) failed, e.g. %s (%s): %s", nrow(failed), failed$gene[1], failed$chr[1], failed$error[1])
    stop(structure(class = c("nma_region_failure", "error", "condition"),
                   list(message = msg, call = NULL, failed = failed)))
  }
  saturated <- rbindlist(lapply(seq_along(res), function(i)
    if (length(res[[i]]$saturated)) data.table(chr = regions$chr[i], pos = res[[i]]$saturated)))

  out <- rbindlist(lapply(res, `[[`, "rows"), use.names = TRUE, fill = TRUE)
  if (nrow(out) == 0) return(list(counts = empty_counts(), failed = failed, duplicate_sites = 0L,
                                  saturated = saturated))
  for (col in BASE_COLUMNS) {
    if (!col %in% names(out)) out[, (col) := 0L]
    out[is.na(get(col)), (col) := 0L]
    out[, (col) := as.integer(get(col))]
  }
  out[, chr := as.character(chr)]
  out[, totalReads := A + C + G + T + Deletion + Insertion]
  out <- out[totalReads > min_coverage]
  out[, delrate := Deletion / totalReads]
  setnames(out, BASE_COLUMNS, paste0(BASE_COLUMNS, ".count"))
  out[, pos := as.integer(pos)]
  out <- out[, ..COUNT_COLUMNS]
  dup <- sum(duplicated(out[, .(chr, pos)]))
  list(counts = out, failed = failed, duplicate_sites = as.integer(dup), saturated = saturated)
}

empty_counts <- function() {
  dt <- data.table(chr = character(), pos = integer(), gene = character())
  for (c in COUNT_COLUMNS[4:10]) dt[, (c) := integer()]
  dt[, `:=`(ref = character(), kmer = character(), strand = character(), delrate = numeric())]
  dt[, ..COUNT_COLUMNS]
}

#' Write a counts table (tab-separated, NA written as "NA").
write_counts <- function(dt, file) {
  fwrite(dt, file, sep = "\t", quote = FALSE, na = "NA")
  invisible(file)
}
