# Tiny reference + BAM fixtures built from SAM text, so every expected count is known exactly.

# 40 nt, + strand amplicon. T at 4, 8, 9, 10, 15, 19, 20, 26, 29, 30, 36, 39, 40 (see t_positions()).
REF_SEQ <- "ACGTACGTTTAGCATGCATTACGGATCCTTAGCAATCGTT"
REF_LEN <- nchar(REF_SEQ)
t_positions <- function(seq = REF_SEQ) which(strsplit(seq, "")[[1]] == "T")

make_fasta <- function(dir, seq = REF_SEQ, name = "amp1") {
  fa <- file.path(dir, "ref.fa")
  writeLines(c(paste0(">", name), seq), fa)
  Rsamtools::indexFa(fa)
  fa
}

# One full-length read starting at position 1, with optional edits:
#   del = position deleted; ins = position after which one base is inserted; sub = list(pos, base)
read_spec <- function(name, flag = 0L, mapq = 60L, del = NULL, ins = NULL, sub = NULL, seq = REF_SEQ) {
  bases <- strsplit(seq, "")[[1]]
  n <- length(bases)
  if (!is.null(sub)) bases[sub[[1]]] <- sub[[2]]
  if (!is.null(del)) {
    cigar <- paste0(if (del > 1) paste0(del - 1, "M"), "1D", if (del < n) paste0(n - del, "M"))
    bases <- bases[-del]
  } else if (!is.null(ins)) {
    cigar <- paste0(ins, "M1I", n - ins, "M")
    bases <- append(bases, "A", after = ins)
  } else {
    cigar <- paste0(n, "M")
  }
  s <- paste(bases, collapse = "")
  data.frame(qname = name, flag = flag, rname = "amp1", pos = 1L, mapq = mapq, cigar = cigar,
             rnext = "*", pnext = 0L, tlen = 0L, seq = s, qual = strrep("I", nchar(s)),
             stringsAsFactors = FALSE)
}

reads <- function(n, prefix, ...) do.call(rbind, lapply(seq_len(n), function(i) read_spec(sprintf("%s%03d", prefix, i), ...)))

make_bam <- function(dir, sam_rows, name = "s1", ref_name = "amp1", ref_len = REF_LEN) {
  sam <- file.path(dir, paste0(name, ".sam"))
  writeLines(c("@HD\tVN:1.6\tSO:unsorted", sprintf("@SQ\tSN:%s\tLN:%d", ref_name, ref_len)), sam)
  utils::write.table(sam_rows, sam, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  bam <- Rsamtools::asBam(sam, file.path(dir, name), overwrite = TRUE, indexDestination = TRUE)
  bam
}

write_bed <- function(dir, rows, name = "targets.bed") {
  bed <- file.path(dir, name)
  utils::write.table(rows, bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  bed
}
bed_row <- function(start, end, strand = "+", chr = "amp1", gene = "AMP1")
  data.frame(chr = chr, start = start, end = end, gene = gene, score = 0, strand = strand)

# Standard design: 30 clean reads, 10 with a deletion at 19, 5 with C at 19, 4 MAPQ-0 reads with a
# deletion at 19 (ignored, min_mapq = 1), 3 reverse-strand reads (ignored on a + region).
standard_reads <- function() {
  rbind(reads(30, "clean"), reads(10, "del", del = 19L), reads(5, "sub", sub = list(19L, "C")),
        reads(4, "mq0", mapq = 0L, del = 19L), reads(3, "rev", flag = 16L, del = 19L))
}
