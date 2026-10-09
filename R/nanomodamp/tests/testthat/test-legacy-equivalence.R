# WP4 testthat case 11: exact equivalence with the legacy site-calling code (plan §7 WP4).
#
# Inputs are the paper's shipped count tables (PUS7regulation2026 @ 4fd285b,
# Figure3/parameter_sweep/BIDdetect_data_invitro_delpos.txt and Figure3/plots/
# BIDdetect_data_incell_delpos.txt). References (fixtures/legacy_equivalence/legacy_rerun/) are
# the outputs of modification_analysis.R and incell_analysis.R rerun unchanged on those inputs
# under R 4.3.2 (see fixtures/legacy_equivalence/README.md); the shipped paper tables are
# checked too (categories).
#
# Tolerances: p-values relative 1e-8 (values below 1e-290 count as equal: readr writes subnormal
# doubles with a wrong exponent, e.g. 6.5187e-315 -> "6.5187e-298"); deltas and dd 1e-12
# absolute; categories, labels, site sets and row order exact.
#
# The in-cellulo part runs on a deterministic subset of sites by default (legacy BH is per site,
# R-30, so every site's result is independent of the others). Set
# NANOMODAMP_FULL_EQUIVALENCE=true to run all sites (CI).

fx <- function(...) testthat::test_path("fixtures", "legacy_equivalence", ...)
read_fx <- function(...) readr::read_tsv(fx(...), show_col_types = FALSE, progress = FALSE)
full_run <- function() identical(tolower(Sys.getenv("NANOMODAMP_FULL_EQUIVALENCE")), "true")

same_num <- function(a, b, rel = 1e-8) {
  both_na <- is.na(a) & is.na(b)
  tiny <- !is.na(a) & !is.na(b) & abs(a) < 1e-290 & abs(b) < 1e-290
  inf <- !is.na(a) & !is.na(b) & is.infinite(a) & is.infinite(b) & a == b
  ok <- !is.na(a) & !is.na(b) & is.finite(a) & is.finite(b) & abs(a - b) <= rel * pmax(abs(b), .Machine$double.xmin)
  both_na | tiny | inf | ok
}
expect_same_num <- function(a, b, what, rel = 1e-8) {
  bad <- which(!same_num(a, b, rel))
  expect(length(bad) == 0, sprintf("%s: %d value(s) differ, e.g. %s vs %s", what, length(bad),
                                   format(a[bad[1]], digits = 15), format(b[bad[1]], digits = 15)))
}
expect_same_abs <- function(a, b, what, tol = 1e-12) {
  bad <- which(!((is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= tol)))
  expect(length(bad) == 0, sprintf("%s: %d value(s) differ", what, length(bad)))
}

# incell_analysis.R uses "Equivalent" / "Not equivalent"; the port uses the detailed labels of
# modification_analysis.R (plan §6.3). Map detailed -> simple for the comparison.
simple_label <- function(x) {
  ifelse(is.na(x) | x %in% c("Not enough data", "Model fitting failed", "Eq. test not run"), x,
         ifelse(x == "Equivalent by TOST", "Equivalent", "Not equivalent"))
}

incell_input <- function() {
  d <- read_fx("inputs", "BIDdetect_data_incell_delpos.txt.gz")
  d$treat <- ifelse(d$treat %in% c("in", "input"), "input", "BS")
  if (!full_run()) {
    # 20 sites the paper called PUS7-dependent (so the positive paths are covered) + the first
    # 40 sites in (chr, pos) order
    called <- dplyr::bind_rows(dplyr::tibble(chr = character(), pos = double()),
                               lapply(Sys.glob(fx("legacy_rerun", "incell", "PUS7_dep_*_significant_summary.tsv")),
                                      function(f) dplyr::select(readr::read_tsv(f, show_col_types = FALSE), chr, pos))) |>
      dplyr::distinct() |> dplyr::arrange(chr, pos) |> utils::head(20)
    others <- dplyr::distinct(d, chr, pos) |> dplyr::arrange(chr, pos) |> utils::head(40)
    d <- dplyr::semi_join(d, dplyr::distinct(dplyr::bind_rows(called, others)), by = c("chr", "pos"))
  }
  d
}

test_that("11a: in vitro treatment test equals modification_analysis.R exactly", {
  skip_if_not(file.exists(fx("inputs", "BIDdetect_data_invitro_delpos.txt.gz")), "fixtures missing")
  d <- read_fx("inputs", "BIDdetect_data_invitro_delpos.txt.gz")
  # D14: the sample sheet labels noPUS samples treat = input (the legacy script recoded them)
  d$treat <- ifelse(d$vector == "noPUS" | d$treat %in% c("in", "input"), "input", "BS")
  r <- suppressWarnings(treatment_test(d))
  ref <- read_fx("legacy_rerun", "invitro", "invitro_delpos_modification_significance.tsv")
  expect_equal(nrow(r), nrow(ref))
  expect_identical(paste(r$chr, r$pos), paste(ref$chr, ref$pos))   # same order (delta desc)
  expect_identical(r$category, ref$category)
  expect_identical(r$equivalence_status, ref$equivalence_status)
  for (col in c("avg_delrate_BS", "avg_delrate_input", "delta_delrate", "site_specific_eqbound")) {
    expect_same_abs(r[[col]], ref[[col]], col)
  }
  expect_same_num(r$p.adjust_diff, ref$p.adjust_diff, "p.adjust_diff")

  shipped <- read_fx("shipped", "invitro_modification_significance.tsv")
  m <- dplyr::inner_join(r, shipped, by = c("chr", "pos"), suffix = c("", ".paper"))
  expect_equal(nrow(m), nrow(shipped))
  expect_identical(m$category, m$category.paper)   # paper categories reproduced
})

wt_mod_runs <- list(
  HepG2 = list(celltype = "HepG2", re = "(1|vector)"),
  `293T` = list(celltype = "293T", re = "(1|vector)"),
  Both = list(celltype = c("HepG2", "293T"), re = "(1|vector) + (1|celltype)"))

test_that("11b: in-cellulo WT_mod treatment tests equal incell_analysis.R", {
  skip_if_not(file.exists(fx("inputs", "BIDdetect_data_incell_delpos.txt.gz")), "fixtures missing")
  d <- incell_input()
  for (cond in names(wt_mod_runs)) {
    run <- wt_mod_runs[[cond]]
    ds <- subset_rows(d, list(celltype = run$celltype, vector = c("WT", "P102", "P4")))
    r <- suppressWarnings(treatment_test(ds, random_effects = run$re))
    ref <- read_fx("legacy_rerun", "incell", paste0("WT_mod_", cond, "_significance.tsv")) |>
      dplyr::semi_join(r, by = c("chr", "pos"))
    m <- dplyr::inner_join(r, ref, by = c("chr", "pos"), suffix = c("", ".ref"))
    expect_equal(nrow(m), nrow(r), label = paste("rows", cond))
    expect_identical(m$category, m$category.ref, label = paste("category", cond))
    expect_identical(simple_label(m$equivalence_status), m$equivalence_status.ref, label = paste("labels", cond))
    expect_identical(m$is_equivalent, m$is_equivalent.ref, label = paste("is_equivalent", cond))
    expect_identical(m$all_below_thresh, m$all_below_thresh.ref, label = paste("prefilter", cond))
    for (col in c("avg_delrate_BS", "avg_delrate_input", "delta_delrate", "site_specific_eqbound")) {
      expect_same_abs(m[[col]], m[[paste0(col, ".ref")]], paste(cond, col))
    }
    expect_same_num(m$p.value, m$p.value.ref, paste(cond, "p.value"))
    expect_same_num(m$p.adjust_diff, m$p.adjust_diff.ref, paste(cond, "p.adjust_diff"))
    if (cond == "Both" && full_run()) {
      # Paper vs today's rerun of the same legacy code: one site changes category because the
      # original fit differed (p 1.4e-5 vs 3.0e-22; package/optimizer drift, not the port).
      paper_drift <- "RPL22_chr1_6186768:Inconclusive>Unmodified"
      shipped <- read_fx("shipped", "WT_mod_Both_significance.tsv")
      mm <- dplyr::inner_join(r, shipped, by = c("chr", "pos"), suffix = c("", ".paper"))
      expect_equal(nrow(mm), nrow(shipped))
      diff <- mm[mm$category != mm$category.paper, ]
      expect_setequal(paste0(diff$chr, ":", diff$category.paper, ">", diff$category), paper_drift)
    }
  }
})

pus7_runs <- list(
  HepG2_WT_v_KD = list(celltype = "HepG2", levels = c("P101", "P102"), re = ""),
  HepG2_WT_v_OE = list(celltype = "HepG2", levels = c("P4", "P3"), re = ""),
  `293T_WT_v_KD` = list(celltype = "293T", levels = c("P101", "P102"), re = ""),
  `293T_WT_v_OE` = list(celltype = "293T", levels = c("P4", "P3"), re = ""),
  Both_WT_v_KD = list(celltype = c("HepG2", "293T"), levels = c("P101", "P102"), re = "(1|celltype)"),
  Both_WT_v_OE = list(celltype = c("HepG2", "293T"), levels = c("P4", "P3"), re = "(1|celltype)"))

test_that("11c: PUS7-dependency factor tests and their union equal incell_analysis.R", {
  skip_if_not(file.exists(fx("inputs", "BIDdetect_data_incell_delpos.txt.gz")), "fixtures missing")
  d <- incell_input()
  sites <- dplyr::distinct(d, chr, pos)
  sig <- list()
  for (run_name in names(pus7_runs)) {
    run <- pus7_runs[[run_name]]
    ds <- subset_rows(d, list(celltype = run$celltype, vector = run$levels))
    f <- suppressWarnings(factor_test(ds, factor = "vector", levels = run$levels, random_effects = run$re))
    sig[[run_name]] <- f$significant
    ref_file <- fx("legacy_rerun", "incell", paste0("PUS7_dep_", run_name, "_significant_summary.tsv"))
    ref <- if (file.exists(ref_file)) readr::read_tsv(ref_file, show_col_types = FALSE) else f$significant[0, ]
    ref <- dplyr::semi_join(ref, sites, by = c("chr", "pos"))
    expect_setequal(paste(f$significant$chr, f$significant$pos), paste(ref$chr, ref$pos))
    m <- dplyr::inner_join(f$significant, ref, by = c("chr", "pos"), suffix = c("", ".ref"))
    expect_identical(m$term, m$term.ref)
    expect_same_num(m$p.value, m$p.value.ref, paste(run_name, "p.value"))
    expect_same_num(m$p.value.BH, m$p.value.BH.ref, paste(run_name, "p.value.BH"))
    for (col in c(paste0("delta_delrate_", run$levels), "dd_delrate")) {
      expect_same_abs(m[[col]], m[[paste0(col, ".ref")]], paste(run_name, col))
    }
  }
  u <- site_sets(sig, "union")
  ref_u <- read_fx("legacy_rerun", "incell", "PUS7_dep_union_significant_summary.tsv") |>
    dplyr::semi_join(sites, by = c("chr", "pos"))
  expect_setequal(paste(u$chr, u$pos), paste(ref_u$chr, ref_u$pos))
  shipped_u <- read_fx("shipped", "PUS7_dep_union_significant_summary.tsv") |> dplyr::semi_join(sites, by = c("chr", "pos"))
  expect_setequal(paste(u$chr, u$pos), paste(shipped_u$chr, shipped_u$pos))  # paper union reproduced
  m <- dplyr::inner_join(u, ref_u, by = c("chr", "pos"), suffix = c("", ".ref"))
  for (run_name in names(pus7_runs)) {
    for (col in paste0(c("p.value_", "dd_delrate_"), run_name)) {
      expect_same_num(m[[col]], m[[paste0(col, ".ref")]], paste("union", col))
    }
  }
})
