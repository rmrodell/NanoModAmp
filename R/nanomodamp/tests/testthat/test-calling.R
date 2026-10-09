# WP4 testthat cases 1–9 (plan §7 WP4) on simulated data; case 10 is in test-plots.R and case 11
# (legacy equivalence) in test-legacy-equivalence.R.

treat_res <- function(...) treatment_test(sim_treatment(), ...)

test_that("1: clear positive is Modified, clear null Unmodified, low-coverage noisy Inconclusive", {
  r <- treat_res()
  cat_of <- function(chr) r$category[r$chr == chr]
  expect_equal(cat_of("POS"), "Modified")
  expect_equal(cat_of("NULL"), "Unmodified")
  expect_equal(r$equivalence_status[r$chr == "NULL"], "Equivalent by TOST")
  expect_equal(cat_of("NOISY"), "Inconclusive")
})

test_that("2: sites with every delrate below SESOI are pre-filtered as Unmodified and not tested", {
  r <- treat_res()
  low <- r[r$chr == "LOW", ]
  expect_true(low$all_below_thresh)
  expect_equal(low$category, "Unmodified")
  expect_true(is.na(low$p.value))
  expect_true(is.na(low$equivalence_status))
})

test_that("3: 'Not enough data' and 'Model fitting failed' paths", {
  d <- sim_treatment()
  one_rep <- d[d$chr == "POS" & d$rep == 1, ]
  r <- treatment_test(one_rep)
  expect_equal(r$equivalence_status, "Not enough data")
  expect_equal(r$category, "Inconclusive")

  local_mocked_bindings(bglmer = function(...) stop("forced failure"), .package = "blme")
  r <- treatment_test(d[d$chr == "POS", ])
  expect_equal(r$equivalence_status, "Model fitting failed")
  expect_true(is.na(r$p.value))
  f <- factor_test(sim_factor(), factor = "level", levels = c("A", "B"))
  expect_true(all(f$all_tests$note == "Model fitting failed"))
  expect_equal(nrow(f$significant), 0)
})

test_that("4: p.adjust gets exactly one p-value per tested site (BH), legacy is per-site (R-30)", {
  d <- sim_treatment()
  bh <- treatment_test(d, p_adjust = "BH")
  tested <- !is.na(bh$p.value)
  expect_equal(sum(tested), 3)  # LOW is pre-filtered
  expect_equal(bh$p.adjust_diff[tested], stats::p.adjust(bh$p.value[tested], method = "BH"))
  legacy <- treatment_test(d)
  expect_equal(legacy$p.adjust_diff, legacy$p.value)

  f <- factor_test(sim_factor(), factor = "level", levels = c("A", "B"), p_adjust = "BH")
  kept <- f$all_tests$term %in% "glm_model_factor"
  expect_equal(sum(kept), 3)
  expect_equal(f$all_tests$p.value.BH[kept], stats::p.adjust(f$all_tests$p.value[kept], method = "BH"))
})

test_that("5: input delrate = 0 behaves as the legacy code (D25: infinite bound)", {
  d <- sim_site("ZERO", 10, c(input = 0, BS = 0.30), n = 500, seed = 21)
  r <- treatment_test(dplyr::as_tibble(d))
  expect_equal(r$avg_delrate_input, 0)
  expect_equal(r$site_specific_eqbound, Inf)
  expect_true(r$is_equivalent)                       # any finite CI lies within (-Inf, Inf)
  expect_equal(r$equivalence_status, "Equivalent by TOST")
  expect_equal(r$category, "Modified")               # Modified is checked before equivalence
})

test_that("6: an interaction is detected for positive and both, and the sign split is correct", {
  d <- sim_factor()
  pos <- factor_test(d, factor = "level", levels = c("A", "B"), direction = "positive")
  expect_equal(pos$significant$chr, "INTER")
  expect_gt(pos$significant$dd_delrate, 0.05)
  expect_named(pos$significant, c("chr", "pos", "term", "p.value", "note", "p.value.BH",
                                  "delta_delrate_A", "delta_delrate_B", "dd_delrate"))
  both <- factor_test(d, factor = "level", levels = c("A", "B"), direction = "both")
  expect_setequal(both$significant$chr, c("INTER", "REVERSE"))
  expect_equal(both$experimental_high$chr, "INTER")
  expect_equal(both$baseline_high$chr, "REVERSE")
})

test_that("7: a site with equal effects in both levels is not significant", {
  f <- factor_test(sim_factor(), factor = "level", levels = c("A", "B"), direction = "both")
  expect_false("EQUAL" %in% f$significant$chr)
  eq <- f$all_tests[f$all_tests$chr == "EQUAL", ]
  expect_lt(abs(eq$dd_delrate), 0.05)
})

test_that("8: seeded runs on 1 and 2 workers are identical", {
  skip_on_cran()
  d <- sim_treatment()
  expect_identical(treatment_test(d, workers = 1), treatment_test(d, workers = 2))
  f <- sim_factor()
  expect_identical(factor_test(f, "level", c("A", "B"), workers = 1),
                   factor_test(f, "level", c("A", "B"), workers = 2))
})

test_that("9: result tables follow the §5.5 contracts", {
  r <- treat_res()
  expect_equal(names(r)[1:9], c("chr", "pos", "category", "avg_delrate_BS", "avg_delrate_input",
                                "delta_delrate", "p.adjust_diff", "equivalence_status",
                                "site_specific_eqbound"))
  expect_true(all(diff(r$delta_delrate) <= 0))  # sorted by delta_delrate descending
  expect_true(all(r$category %in% c("Modified", "Unmodified", "Inconclusive")))
  f <- factor_test(sim_factor(), factor = "level", levels = c("A", "B"))
  cols <- c("chr", "pos", "term", "p.value", "note", "p.value.BH", "delta_delrate_A",
            "delta_delrate_B", "dd_delrate")
  expect_named(f$all_tests, cols)
  expect_named(f$significant, cols)
  expect_true(all(f$significant$term == "glm_model_factor"))
})

test_that("input validation: explicit treat values (D14) and required columns", {
  d <- sim_treatment()
  d$treat[1] <- "noPUS"
  expect_error(treatment_test(d), "treat must be 'input' or 'BS'")
  expect_error(treatment_test(d[, setdiff(names(d), "rep")]), "missing required columns: rep")
  expect_error(factor_test(sim_factor(), "nope", c("A", "B")), "not in the counts table")
})

test_that("site sets: union, intersection, difference and quartiles", {
  a <- dplyr::tibble(chr = c("x", "y"), pos = c(1, 2), dd_delrate = c(0.1, 0.2))
  b <- dplyr::tibble(chr = c("y", "z"), pos = c(2, 3), dd_delrate = c(0.3, 0.4))
  u <- site_sets(list(a = a, b = b), "union")
  expect_equal(nrow(u), 3)
  expect_named(u, c("chr", "pos", "dd_delrate_a", "dd_delrate_b"))
  expect_equal(nrow(site_sets(list(a = a, b = b), "intersection")), 1)
  expect_equal(site_sets(list(a = a, b = b), "difference")$chr, "x")
  ref <- dplyr::tibble(chr = letters[1:8], pos = 1, delta_delrate = 1:8 / 10)
  q <- site_quartiles(dplyr::tibble(chr = letters[1:8], pos = 1), ref, "delta_delrate")
  expect_equal(q$quartile, rep(1:4, each = 2))
  expect_null(site_quartiles(dplyr::tibble(chr = letters[1:3], pos = 1), ref, "delta_delrate"))
})
