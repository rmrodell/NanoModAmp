# Legacy equivalence: count_sites() vs the paper's counter on the same fixture BAM.
#
# inst/legacy/bam_counts_fast.R is a byte-identical copy of PUS7regulation2026 @ 4fd285b
# Figure2/BIDdetect/bam_counts_fast.R (md5 d5e1ae8f21fd90c05a8cf7e5825b453b). The legacy code reads
# BED start as-is, so count_sites() is run with bed_coordinates = "one_based_start". Expected,
# documented differences: column `reference` is `ref` (§5.3) and kmer NA is written as "NA"
# instead of an empty field. Golden-package BAMs are not shipped, hence the fixture BAM.

legacy_deps <- c("future", "future.apply", "future.callr", "optparse", "data.table", "GenomicRanges",
                 "GenomicAlignments", "Rsamtools", "Biostrings", "BiocManager")

test_that("count_sites() equals the legacy bam_counts_fast.R on the same BAM", {
  missing <- legacy_deps[!vapply(legacy_deps, requireNamespace, logical(1), quietly = TRUE)]
  skip_if(length(missing) > 0, paste("legacy counter needs", paste(missing, collapse = ", ")))
  legacy <- system.file("legacy", "bam_counts_fast.R", package = "nanomodamp")
  skip_if(legacy == "", "legacy script not installed")

  d <- tempfile("nmaleg"); dir.create(d)
  fa <- make_fasta(d)
  rows <- rbind(standard_reads(), reads(3, "ins", ins = 9L), reads(2, "sub8", sub = list(8L, "G")),
                reads(4, "rcd", flag = 16L, del = 22L), reads(25, "rc", flag = 16L))
  bam <- make_bam(d, rows)
  bed <- write_bed(d, rbind(bed_row(1, REF_LEN, gene = "FULL"), bed_row(19, 19, gene = "SITE"),
                            bed_row(1, REF_LEN, strand = "-", gene = "MINUS")))
  out <- file.path(d, "legacy.tsv")
  rscript <- file.path(R.home("bin"), "Rscript")
  st <- system2(rscript, c(legacy, "--bedFile", bed, "--bamFile", bam, "--referenceFasta", fa,
                           "--outputFile", out, "--threads", "1"), stdout = FALSE, stderr = FALSE)
  expect_identical(st, 0L)
  leg <- data.table::fread(out, sep = "\t", na.strings = c("", "NA"))
  data.table::setnames(leg, "reference", "ref")

  new <- count_sites(bam, fa, bed, bed_coordinates = "one_based_start")$counts
  expect_identical(names(leg), names(new))
  expect_gt(nrow(new), 20)
  for (col in names(new)) {
    expect_equal(as.vector(leg[[col]]), as.vector(new[[col]]), info = col)
  }
})
