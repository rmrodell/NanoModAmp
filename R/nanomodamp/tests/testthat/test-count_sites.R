# WP3 testthat cases 1-9, 11, 12 (plan §7 WP3) plus R-12 and R-16.

fixture <- function(rows = standard_reads(), seq = REF_SEQ) {
  d <- withr_tempdir()
  list(dir = d, fa = make_fasta(d, seq), bam = make_bam(d, rows))
}
withr_tempdir <- function() { d <- tempfile("nma"); dir.create(d); d }
site <- function(res, p) res$counts[pos == p]

test_that("1: deletion count equals the truth (MAPQ 0 and wrong-strand reads ignored)", {
  f <- fixture()
  res <- count_sites(f$bam, f$fa, write_bed(f$dir, bed_row(0, REF_LEN)))
  s <- site(res, 19)
  expect_equal(nrow(s), 1)
  expect_identical(s$Deletion.count, 10L)
  expect_identical(s$T.count, 30L)
  expect_identical(s$C.count, 5L)
  expect_identical(s$totalReads, 45L)
  expect_equal(s$delrate, 10 / 45)
  expect_identical(names(res$counts), COUNT_COLUMNS)
})

test_that("2: only T reference positions are reported by default", {
  f <- fixture()
  res <- count_sites(f$bam, f$fa, write_bed(f$dir, bed_row(0, REF_LEN)))
  expect_true(all(res$counts$ref == "T"))
  expect_setequal(res$counts$pos, t_positions())
})

test_that("3: all-bases mode reports every covered position", {
  f <- fixture()
  res <- count_sites(f$bam, f$fa, write_bed(f$dir, bed_row(0, REF_LEN)), count_all_bases = TRUE)
  expect_setequal(res$counts$pos, seq_len(REF_LEN))
  expect_identical(res$counts[pos == 1]$ref, "A")
})

test_that("4: kmer is NA within 2 nt of either end and the centred 5-mer elsewhere", {
  f <- fixture()
  res <- count_sites(f$bam, f$fa, write_bed(f$dir, bed_row(0, REF_LEN)), count_all_bases = TRUE)
  expect_true(all(is.na(res$counts[pos %in% c(1, 2, REF_LEN - 1, REF_LEN)]$kmer)))
  expect_identical(res$counts[pos == 19]$kmer, substr(REF_SEQ, 17, 21))
  expect_true(all(!is.na(res$counts[pos %in% 3:(REF_LEN - 2)]$kmer)))
})

test_that("5: a - strand region is complemented and reverse-complemented", {
  rows <- rbind(reads(25, "rc", flag = 16L), reads(5, "rcdel", flag = 16L, del = 19L))
  f <- fixture(rows)
  res <- count_sites(f$bam, f$fa, write_bed(f$dir, bed_row(0, REF_LEN, strand = "-")))
  # ref is complemented: genomic A becomes "T" on the - strand
  a_pos <- which(strsplit(REF_SEQ, "")[[1]] == "A")
  expect_setequal(res$counts$pos, a_pos)
  expect_true(all(res$counts$strand == "-"))
  p <- a_pos[a_pos >= 3 & a_pos <= REF_LEN - 2][1]
  expect_identical(res$counts[pos == p]$kmer,
                   as.character(Biostrings::reverseComplement(Biostrings::DNAString(substr(REF_SEQ, p - 2, p + 2)))))
  expect_identical(res$counts[pos == p]$T.count, 30L)  # genomic A read on - strand counted as T
  # the deletion at genomic 19 (a T, i.e. an A on -) is not a reported site
  expect_false(19 %in% res$counts$pos)
})

test_that("6: reads on the other strand are excluded", {
  f <- fixture(rbind(reads(25, "fwd"), reads(40, "rev", flag = 16L, del = 19L)))
  res <- count_sites(f$bam, f$fa, write_bed(f$dir, bed_row(0, REF_LEN)))
  expect_identical(site(res, 19)$Deletion.count, 0L)
  expect_identical(site(res, 19)$totalReads, 25L)
})

test_that("7: single-site and amplicon BEDs under both bed_coordinates conventions", {
  f <- fixture()
  one <- function(rows, conv) count_sites(f$bam, f$fa, write_bed(f$dir, rows, "x.bed"), bed_coordinates = conv)$counts$pos
  expect_identical(one(bed_row(18, 19), "bed0"), 19L)               # standard single site
  expect_identical(one(bed_row(19, 19), "one_based_start"), 19L)    # legacy single site (paper delpos BEDs)
  expect_identical(one(bed_row(18, 19), "one_based_start"), c(19L))  # 18 is not a T: no extra row
  expect_identical(one(bed_row(19, 20), "one_based_start"), c(19L, 20L))  # legacy reads start as 1-based
  expect_identical(one(bed_row(19, 20), "bed0"), 20L)
  expect_setequal(one(bed_row(0, REF_LEN), "bed0"), t_positions())
  expect_setequal(one(bed_row(1, REF_LEN), "one_based_start"), t_positions())
  expect_error(one(bed_row(19, 19), "bed0"), "start = end")
})

test_that("8: an insertion next to a T adds exactly 1 to that position's totalReads (D9)", {
  f0 <- fixture(reads(31, "c"))                                   # 31 reads without an insertion
  f1 <- fixture(rbind(reads(30, "c"), reads(1, "ins", ins = 19L)))  # 30 + 1 read with an insertion after 19
  b0 <- write_bed(f0$dir, bed_row(0, REF_LEN)); b1 <- write_bed(f1$dir, bed_row(0, REF_LEN))
  r0 <- count_sites(f0$bam, f0$fa, b0)$counts; r1 <- count_sites(f1$bam, f1$fa, b1)$counts
  m <- merge(r0, r1, by = "pos", suffixes = c(".0", ".1"))
  d <- m[totalReads.1 != totalReads.0]
  expect_equal(nrow(d), 1)
  expect_identical(d$totalReads.1 - d$totalReads.0, 1L)
  expect_identical(d$Insertion.count.1, 1L)
  expect_identical(d$pos, 19L)  # Rsamtools reports the insertion at the base before it
})

test_that("9: a site with exactly 20 reads is excluded and one with 21 is kept (D10)", {
  f20 <- fixture(reads(20, "a"))
  f21 <- fixture(reads(21, "a"))
  expect_equal(nrow(count_sites(f20$bam, f20$fa, write_bed(f20$dir, bed_row(18, 19)))$counts), 0)
  expect_equal(nrow(count_sites(f21$bam, f21$fa, write_bed(f21$dir, bed_row(18, 19)))$counts), 1)
})

test_that("11: rerunning gives identical output files", {
  f <- fixture()
  bed <- write_bed(f$dir, bed_row(0, REF_LEN))
  o1 <- file.path(f$dir, "a.tsv"); o2 <- file.path(f$dir, "b.tsv")
  write_counts(count_sites(f$bam, f$fa, bed)$counts, o1)
  write_counts(count_sites(f$bam, f$fa, bed)$counts, o2)
  expect_identical(unname(tools::md5sum(o1)), unname(tools::md5sum(o2)))
})

test_that("12: 1 thread and N threads give identical output", {
  f <- fixture()
  bed <- write_bed(f$dir, rbind(bed_row(0, 20, gene = "A"), bed_row(20, REF_LEN, gene = "B"),
                                bed_row(18, 19, gene = "C")))
  expect_identical(count_sites(f$bam, f$fa, bed, threads = 1)$counts,
                   count_sites(f$bam, f$fa, bed, threads = 2)$counts)
})

test_that("overlapping regions repeat rows as in legacy and are counted (§6.2 step 8)", {
  f <- fixture()
  res <- count_sites(f$bam, f$fa, write_bed(f$dir, rbind(bed_row(0, REF_LEN), bed_row(18, 19, gene = "S"))))
  expect_identical(res$duplicate_sites, 1L)
})

test_that("R-12: a failing region fails the task unless allowed, and is listed", {
  f <- fixture()
  bed <- write_bed(f$dir, rbind(bed_row(0, REF_LEN), bed_row(0, 10, chr = "missing", gene = "M")))
  expect_error(count_sites(f$bam, f$fa, bed), class = "nma_region_failure")
  res <- count_sites(f$bam, f$fa, bed, allow_region_failures = TRUE)
  expect_identical(res$failed$gene, "M")
  expect_true(nrow(res$counts) > 0)
})

test_that("R-16: reaching max_depth is reported", {
  f <- fixture()
  res <- count_sites(f$bam, f$fa, write_bed(f$dir, bed_row(18, 19)), max_depth = 10L, min_coverage = 0)
  expect_true(19L %in% res$saturated$pos)
  expect_equal(nrow(count_sites(f$bam, f$fa, write_bed(f$dir, bed_row(18, 19)))$saturated), 0)
})

test_that("R-31: a stale .fai (same names and lengths, wrong offsets) is rejected; a fresh one is used", {
  d <- tempfile("nmastale"); dir.create(d)
  fa <- file.path(d, "ref.fa")
  # index a CRLF version, then rewrite the FASTA with LF line endings, keeping the old index
  writeBin(charToRaw(paste0(">amp1\r\n", REF_SEQ, "\r\n")), fa)
  Rsamtools::indexFa(fa)
  writeLines(c(">amp1", REF_SEQ), fa)
  stale <- read.delim(paste0(fa, ".fai"), header = FALSE)
  expect_identical(stale$V2, REF_LEN)                      # names and lengths still match
  bam <- make_bam(d, standard_reads())
  bed <- write_bed(d, bed_row(18, 19))
  expect_error(count_sites(bam, fa, bed), "stale")
  Rsamtools::indexFa(fa)                                   # what SAMTOOLS_FAIDX does in the pipeline
  expect_identical(count_sites(bam, fa, bed)$counts$Deletion.count, 10L)
})
