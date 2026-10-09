# WP4 testthat case 10 (plots: every expected file exists, non-empty, as PDF and PNG) and the
# analyses-YAML / runner layer used by bin/nma_call.R and bin/nma_sitesets.R.

write_cfg <- function(text) {
  f <- withr::local_tempfile(fileext = ".yaml", .local_envir = parent.frame())
  writeLines(text, f)
  f
}

pdf_png <- function(dir, stem) file.path(dir, paste0(stem, c(".pdf", ".png")))

test_that("10: every expected plot exists, is non-empty and comes as PDF and PNG", {
  out <- withr::local_tempdir()
  d <- sim_treatment()
  r <- treatment_test(d)
  # duplicate the Modified site so the quartile plot (needs >= 4 Modified sites) is drawn
  four <- dplyr::bind_rows(lapply(1:4, function(i) dplyr::mutate(r[r$category == "Modified", ], chr = paste0("M", i))))
  files <- plot_treatment(dplyr::bind_rows(r, four), d, out, "t", plot_all_sites = TRUE)
  stems <- c("t_boxplot_modified_summary", "t_heatmap_modified_summary_nolabels",
             "t_heatmap_modified_summary_labels", "t_boxplot_modified_quartiles", "t_allsites")
  expected <- unlist(lapply(stems, function(s) pdf_png(out, s)))
  expect_setequal(files, expected)
  expect_true(all(file.exists(expected)))
  expect_true(all(file.size(expected) > 0))

  f <- factor_test(sim_factor(), factor = "level", levels = c("A", "B"))
  files <- plot_factor(f$significant, sim_factor(), "level", c("A", "B"), out, "f", plot_all_sites = TRUE)
  expected <- unlist(lapply(c("f_summary_boxplot", "f_summary_heatmap_nolabels", "f_summary_heatmap_labels",
                              "f_allsites"), function(s) pdf_png(out, s)))
  expect_setequal(files, expected)
  expect_true(all(file.size(expected) > 0))
})

test_that("read_analyses validates the §5.4 config, including quartiles_by", {
  ok <- write_cfg(c(
    "analyses:",
    "  - {name: t, type: treatment}",
    "  - {name: f, type: factor, factor: level, levels: [A, B]}",
    "  - {name: g, type: factor, factor: level, levels: [A, B], direction: both}",
    "site_sets:",
    "  - {name: u, op: union, of: [f, g], quartiles_by: {analysis: t, column: delta_delrate}}"))
  cfg <- read_analyses(ok)
  expect_named(cfg$analyses, c("t", "f", "g"))
  expect_equal(cfg$analyses$t$sesoi, 0.05)
  expect_equal(cfg$analyses$f$direction, "positive")

  bad_q <- write_cfg(c("analyses:", "  - {name: t, type: treatment}", "  - {name: f, type: factor, factor: level, levels: [A, B]}",
                       "  - {name: g, type: factor, factor: level, levels: [A, B]}",
                       "site_sets:", "  - {name: u, op: union, of: [f, g], quartiles_by: {analysis: nope, column: delta_delrate}}"))
  expect_error(read_analyses(bad_q), "quartiles_by.analysis 'nope' is not a defined analysis")
  bad_of <- write_cfg(c("analyses:", "  - {name: t, type: treatment}", "site_sets:", "  - {name: u, op: union, of: [t, x]}"))
  expect_error(read_analyses(bad_of), "unknown analyses: x")
  bad_key <- write_cfg(c("analyses:", "  - {name: t, type: treatment, sesio: 0.1}"))
  expect_error(read_analyses(bad_key), "unknown keys sesio")
  bad_factor <- write_cfg(c("analyses:", "  - {name: f, type: factor, factor: level}"))
  expect_error(read_analyses(bad_factor), "needs 'factor' and two 'levels'")
})

test_that("the shipped analyses configs are valid", {
  root <- testthat::test_path("..", "..", "..", "..")
  files <- c(file.path(root, "assets", "analyses_example.yaml"), Sys.glob(file.path(root, "assets", "analyses", "*.yaml")))
  skip_if(!file.exists(files[1]), "repo assets not available")
  for (f in files) expect_type(read_analyses(f)$analyses, "list")
})

test_that("run_analysis and run_site_sets write the §4 layout", {
  out <- withr::local_tempdir()
  d <- bind_sites(sim_treatment() |> dplyr::mutate(level = "A"), sim_factor())
  cfg_file <- write_cfg(c(
    "analyses:",
    "  - {name: trt, type: treatment, subset: {chr: [POS, NULL, LOW, NOISY]}}",
    "  - {name: fac, type: factor, subset: {chr: [INTER, EQUAL, REVERSE]}, factor: level, levels: [A, B]}",
    "  - {name: fac2, type: factor, subset: {chr: [INTER, EQUAL, REVERSE]}, factor: level, levels: [A, B], direction: both}",
    "site_sets:",
    "  - {name: both_sets, op: union, of: [fac, fac2]}"))
  cfg <- read_analyses(cfg_file)
  for (a in cfg$analyses) run_analysis(d, a, file.path(out, "calling"))
  expect_true(file.exists(file.path(out, "calling", "trt", "data_summary", "trt_modification_significance.tsv")))
  expect_true(file.exists(file.path(out, "calling", "trt", "trt_log.txt")))
  expect_true(file.exists(file.path(out, "calling", "trt", "trt_resolved_config.yaml")))
  expect_true(file.exists(file.path(out, "calling", "fac", "data_summary", "fac_all_tests.tsv")))
  expect_true(file.exists(file.path(out, "calling", "fac2", "data_summary", "fac2_A_high_summary.tsv")))
  run_site_sets(cfg, file.path(out, "calling"), file.path(out, "site_sets"), d)
  u <- readr::read_tsv(file.path(out, "site_sets", "both_sets", "data_summary", "both_sets_summary.tsv"), show_col_types = FALSE)
  expect_setequal(u$chr, c("INTER", "REVERSE"))
  expect_error(run_analysis(d, list(name = "empty", type = "treatment", subset = list(chr = "none"), random_effects = "",
                                    sesoi = 0.05, fdr = 0.05, plot_all_sites = FALSE), out),
               "subset selects no rows")
})
