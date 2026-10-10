# WP3 case 10 (merge by name) and D32 (--input_counts) cases.

tmp <- function() { d <- tempfile("nmamerge"); dir.create(d); d }

counts_table <- function(n = 2, shift = 0L) {
  data.table::data.table(chr = "amp1", pos = 19L + shift + seq_len(n) - 1L, gene = "AMP1",
                         totalReads = 45L, A.count = 0L, C.count = 5L, G.count = 0L, T.count = 30L,
                         Deletion.count = 10L, Insertion.count = 0L, ref = "T", kmer = NA_character_,
                         strand = "+", delrate = 10 / 45)
}
write_ss <- function(d, ids, extra = list()) {
  ss <- data.table::data.table(sample_id = ids, fastq = paste0("fastq/", ids, ".fastq.gz"),
                               treat = ifelse(grepl("BS", ids), "BS", "input"), rep = "1")
  for (n in names(extra)) ss[, (n) := extra[[n]]]
  f <- file.path(d, "samplesheet.csv"); data.table::fwrite(ss, f); f
}

test_that("10: merging per-sample tables with shuffled column order gives the same result", {
  d <- tmp()
  ss <- write_ss(d, c("s_BS", "s_input"), list(celltype = c("HepG2", "HepG2")))
  a <- file.path(d, "s_BS.counts.tsv"); b <- file.path(d, "s_input.counts.tsv")
  write_counts(counts_table(), a); write_counts(counts_table(), b)
  m1 <- merge_counts(c(a, b), ss)
  t <- counts_table(); write_counts(t[, rev(names(t)), with = FALSE], b)  # shuffled columns
  m2 <- merge_counts(c(b, a), ss)                                         # and file order
  expect_identical(m1, m2)
  expect_identical(names(m1), c("sample_id", "treat", "rep", "celltype", COUNT_COLUMNS))
  expect_identical(unique(m1$sample_id), c("s_BS", "s_input"))           # sample-sheet order
})

test_that("merge_counts rejects count files for unknown samples", {
  d <- tmp()
  ss <- write_ss(d, "a_BS")
  f <- file.path(d, "zz.counts.tsv"); write_counts(counts_table(), f)
  expect_error(merge_counts(f, ss), "not in the sample sheet")
})

merged_table <- function(d, ids, name, meta = list(treat = "BS")) {
  rows <- data.table::rbindlist(lapply(ids, function(i) cbind(data.table::data.table(sample_id = i), as.data.frame(meta), counts_table())))
  f <- file.path(d, name); write_counts(rows, f); f
}

test_that("D32: two tables merge into the same result as one table with all samples", {
  d <- tmp()
  t1 <- merged_table(d, c("a", "b"), "t1.tsv"); t2 <- merged_table(d, "c", "t2.tsv")
  all3 <- merged_table(d, c("a", "b", "c"), "all.tsv")
  res <- merge_count_tables(c(t1, t2), c("run1", "run2"))
  expect_identical(res$merged, nanomodamp:::read_text_table(all3))
  expect_identical(res$sources$source_table, c("run1", "run1", "run2"))
})

test_that("D32: shuffled column order across tables gives the same merge", {
  d <- tmp()
  t1 <- merged_table(d, c("a", "b"), "t1.tsv"); t2 <- merged_table(d, "c", "t2.tsv")
  x <- nanomodamp:::read_text_table(t2)
  t2s <- file.path(d, "t2s.tsv")
  write_counts(x[, c("sample_id", "treat", rev(COUNT_COLUMNS)), with = FALSE], t2s)
  expect_identical(merge_count_tables(c(t1, t2))$merged, merge_count_tables(c(t1, t2s))$merged)
})

test_that("D32: duplicate sample_id across tables fails", {
  d <- tmp()
  expect_error(merge_count_tables(c(merged_table(d, "a", "t1.tsv"), merged_table(d, "a", "t2.tsv"))),
               "more than one table")
})

test_that("D32: differing count/site columns fail", {
  d <- tmp()
  t1 <- merged_table(d, "a", "t1.tsv")
  x <- nanomodamp:::read_text_table(merged_table(d, "b", "t2.tsv"))
  x[, kmer := NULL]
  t2 <- file.path(d, "bad.tsv"); write_counts(x, t2)
  expect_error(merge_count_tables(c(t1, t2)), "count/site columns differ")
})

test_that("D32: metadata columns are the union in first-seen order; missing values become NA with a warning", {
  d <- tmp()
  t1 <- merged_table(d, "a", "t1.tsv", list(treat = "BS", celltype = "HepG2"))
  t2 <- merged_table(d, "b", "t2.tsv", list(treat = "input", run = "r2"))
  w <- capture_warnings(res <- merge_count_tables(c(t1, t2)))
  expect_length(w, 2)
  expect_true(all(grepl("filled with NA", w)))
  expect_identical(names(res$merged)[1:4], c("sample_id", "treat", "celltype", "run"))
  expect_true(all(is.na(res$merged[sample_id == "a"]$run)))
  expect_true(all(is.na(res$merged[sample_id == "b"]$celltype)))
})

test_that("background_qc writes per-sample and per-treatment summaries and both plot formats", {
  d <- tmp()
  m <- data.table::rbindlist(list(
    data.table::data.table(sample_id = "a", treat = "BS", ref = c("A", "C", "T"), delrate = c(0.01, 0.02, 0.3)),
    data.table::data.table(sample_id = "b", treat = "input", ref = c("A", "C", "T"), delrate = c(0.01, 0.01, 0.02))))
  out <- file.path(d, "background")
  res <- background_qc(m, out)
  expect_true(all(file.exists(file.path(out, c("background_qc.tsv", "background_qc_by_treat.tsv",
                                               "background_qc.pdf", "background_qc.png")))))
  expect_equal(nrow(res$per_sample), 6)
  expect_setequal(res$by_treat$treat, c("BS", "input"))
})
