# WP4 testthat case 11: exact equivalence with the legacy site-calling code (plan §7 WP4).
#
# Like with like: at test time the legacy scripts (fixtures/legacy_equivalence/legacy_scripts/,
# byte-identical copies of PUS7regulation2026 @ 4fd285b Figure3/mpra_sites/, md5s in the README)
# and the port run on the same input, in the same R installation and library. Results of GLMM
# fits (bglmer) depend on R / lme4 / blme versions, so a comparison against outputs produced in
# another environment can't hold to 1e-8; a same-environment comparison can.
#
# Inputs are the paper's shipped count tables (Figure3/parameter_sweep/BIDdetect_data_invitro_delpos.txt
# and Figure3/plots/BIDdetect_data_incell_delpos.txt). All port runs use p_adjust = "legacy" (the
# paper's per-site no-op, R-30); the package default is BH.
#
# Tolerances: p-values relative 1e-8 (values below 1e-290 count as equal: readr writes subnormal
# doubles with a wrong exponent, e.g. 6.5187e-315 -> "6.5187e-298"); deltas and dd 1e-12
# absolute; categories, labels, site sets and row order exact.
#
# In cellulo runs on a deterministic subset of sites by default (incell_subset_chrs(), ~50 oligos);
# NANOMODAMP_FULL_EQUIVALENCE=true runs all 760 (legacy ~50 min + port ~35 min on one core).
#
# The stored Sherlock rerun (legacy_rerun/, R 4.3.2, lme4 1.1.35.5, blme 1.0.6, car 3.1.3) and the
# shipped paper tables (shipped/) remain as references: they are asserted only when the running
# environment has those versions, otherwise differences are reported with message().

fx <- function(...) testthat::test_path("fixtures", "legacy_equivalence", ...)
read_fx <- function(...) readr::read_tsv(fx(...), show_col_types = FALSE, progress = FALSE)
read_out <- function(f) readr::read_tsv(f, show_col_types = FALSE, progress = FALSE)
full_run <- function() identical(tolower(Sys.getenv("NANOMODAMP_FULL_EQUIVALENCE")), "true")

# ---- environment ---------------------------------------------------------------------------

# Packages the legacy scripts load that are not nanomodamp imports (the scripts would otherwise try
# to install them).
legacy_pkgs <- c("argparse", "dplyr", "readr", "car", "tidyr", "purrr", "ggplot2", "blme",
                 "future", "furrr", "forcats")
skip_without_legacy_pkgs <- function() {
  for (p in legacy_pkgs) {
    skip_if_not(requireNamespace(p, quietly = TRUE),
                sprintf("the legacy scripts need R package '%s' (not a nanomodamp import)", p))
  }
}

# Versions the stored Sherlock rerun was produced with (fixtures README).
reference_env <- c(lme4 = "1.1.35.5", blme = "1.0.6", car = "3.1.3")
matches_reference_env <- function() {
  getRversion() == "4.3.2" &&
    all(vapply(names(reference_env), function(p) {
      requireNamespace(p, quietly = TRUE) && utils::packageVersion(p) == reference_env[[p]]
    }, logical(1)))
}

# Run a legacy script with the Rscript and library of this session.
run_legacy <- function(script, args = character(), wd) {
  force(script)  # evaluate (relative) paths before changing directory
  force(args)
  libs <- paste(.libPaths(), collapse = .Platform$path.sep)
  old <- setwd(wd)
  on.exit(setwd(old), add = TRUE)
  out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), c(shQuote(script), args),
                                  stdout = TRUE, stderr = TRUE,
                                  env = c(paste0("R_LIBS=", shQuote(libs)), "R_LIBS_USER=NULL")))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    stop("legacy script failed (exit ", status, "):\n", paste(utils::tail(out, 30), collapse = "\n"))
  }
  invisible(out)
}

# incell_analysis.R is configured by editing its USER CONFIGURATION block. Edit a copy at test time
# (the vendored file stays byte-identical) and check that every edit hit exactly one line.
configure_incell <- function(input, output, dest) {
  lines <- sub("\r$", "", readLines(fx("legacy_scripts", "incell_analysis.R"), warn = FALSE))
  edits <- list(
    "^input_file\\s*<-"             = sprintf("input_file <- %s", deparse(input)),
    "^output_dir\\s*<-"             = sprintf("output_dir <- %s", deparse(output)),
    "^num_cores\\s*<-"              = "num_cores <- 1",
    "^plot_all_sites\\s*<-"         = "plot_all_sites <- FALSE",
    "^run_cell_type_analysis\\s*<-" = "run_cell_type_analysis <- FALSE")
  for (pat in names(edits)) {
    hit <- grep(pat, lines)
    stopifnot(length(hit) == 1)
    lines[hit] <- edits[[pat]]
  }
  writeLines(lines, dest)
  dest
}

# ---- comparison helpers (return problem descriptions) -------------------------------------

same_num <- function(a, b, rel = 1e-8) {
  both_na <- is.na(a) & is.na(b)
  tiny <- !is.na(a) & !is.na(b) & abs(a) < 1e-290 & abs(b) < 1e-290
  inf <- !is.na(a) & !is.na(b) & is.infinite(a) & is.infinite(b) & a == b
  ok <- !is.na(a) & !is.na(b) & is.finite(a) & is.finite(b) & abs(a - b) <= rel * pmax(abs(b), .Machine$double.xmin)
  both_na | tiny | inf | ok
}
num_problem <- function(a, b, what, rel = 1e-8) {
  bad <- which(!same_num(a, b, rel))
  if (length(bad)) sprintf("%s: %d value(s) differ, e.g. %s vs %s", what, length(bad),
                           format(a[bad[1]], digits = 15), format(b[bad[1]], digits = 15))
}
abs_problem <- function(a, b, what, tol = 1e-12) {
  bad <- which(!((is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= tol)))
  if (length(bad)) sprintf("%s: %d value(s) differ", what, length(bad))
}
same_problem <- function(a, b, what) {
  if (!identical(a, b)) {
    n <- if (length(a) == length(b)) sum(a != b | xor(is.na(a), is.na(b)), na.rm = TRUE) else NA
    sprintf("%s: %s of %d differ", what, n, length(b))
  }
}

# incell_analysis.R uses "Equivalent" / "Not equivalent"; the port uses the detailed labels of
# modification_analysis.R (plan §6.3). Map detailed -> simple for the comparison.
simple_label <- function(x) {
  ifelse(is.na(x) | x %in% c("Not enough data", "Model fitting failed", "Eq. test not run"), x,
         ifelse(x == "Equivalent by TOST", "Equivalent", "Not equivalent"))
}

diff_invitro <- function(r, ref) {
  if (nrow(r) != nrow(ref)) return(sprintf("rows: %d vs %d", nrow(r), nrow(ref)))
  c(same_problem(paste(r$chr, r$pos), paste(ref$chr, ref$pos), "site order"),
    same_problem(r$category, ref$category, "category"),
    same_problem(r$equivalence_status, ref$equivalence_status, "equivalence_status"),
    unlist(lapply(c("avg_delrate_BS", "avg_delrate_input", "delta_delrate", "site_specific_eqbound"),
                  function(col) abs_problem(r[[col]], ref[[col]], col))),
    num_problem(r$p.adjust_diff, ref$p.adjust_diff, "p.adjust_diff"))
}

diff_wt_mod <- function(r, ref, cond) {
  m <- dplyr::inner_join(r, ref, by = c("chr", "pos"), suffix = c("", ".ref"))
  c(if (nrow(m) != nrow(r) || nrow(m) != nrow(ref))
      sprintf("%s rows: port %d, legacy %d, shared %d", cond, nrow(r), nrow(ref), nrow(m)),
    same_problem(m$category, m$category.ref, paste(cond, "category")),
    same_problem(simple_label(m$equivalence_status), m$equivalence_status.ref, paste(cond, "labels")),
    same_problem(m$is_equivalent, m$is_equivalent.ref, paste(cond, "is_equivalent")),
    same_problem(m$all_below_thresh, m$all_below_thresh.ref, paste(cond, "prefilter")),
    unlist(lapply(c("avg_delrate_BS", "avg_delrate_input", "delta_delrate", "site_specific_eqbound"),
                  function(col) abs_problem(m[[col]], m[[paste0(col, ".ref")]], paste(cond, col)))),
    num_problem(m$p.value, m$p.value.ref, paste(cond, "p.value")),
    num_problem(m$p.adjust_diff, m$p.adjust_diff.ref, paste(cond, "p.adjust_diff")))
}

diff_factor <- function(sig, ref, run_name, levels) {
  m <- dplyr::inner_join(sig, ref, by = c("chr", "pos"), suffix = c("", ".ref"))
  c(same_problem(sort(paste(sig$chr, sig$pos)), sort(paste(ref$chr, ref$pos)), paste(run_name, "significant sites")),
    same_problem(m$term, m$term.ref, paste(run_name, "term")),
    num_problem(m$p.value, m$p.value.ref, paste(run_name, "p.value")),
    num_problem(m$p.value.BH, m$p.value.BH.ref, paste(run_name, "p.value.BH")),
    unlist(lapply(c(paste0("delta_delrate_", levels), "dd_delrate"),
                  function(col) abs_problem(m[[col]], m[[paste0(col, ".ref")]], paste(run_name, col)))))
}

diff_union <- function(u, ref_u, run_names) {
  m <- dplyr::inner_join(u, ref_u, by = c("chr", "pos"), suffix = c("", ".ref"))
  c(same_problem(sort(paste(u$chr, u$pos)), sort(paste(ref_u$chr, ref_u$pos)), "union sites"),
    unlist(lapply(run_names, function(rn) c(
      num_problem(m[[paste0("p.value_", rn)]], m[[paste0("p.value_", rn, ".ref")]], paste("union p.value", rn)),
      abs_problem(m[[paste0("dd_delrate_", rn)]], m[[paste0("dd_delrate_", rn, ".ref")]], paste("union dd", rn))))))
}

# Assert problems against the same-environment legacy run; for the stored references assert only
# in the reference environment, otherwise report.
expect_no_problems <- function(problems, what) {
  expect(length(problems) == 0, paste0(what, ":\n  ", paste(problems, collapse = "\n  ")))
}
check_stored <- function(problems, what) {
  if (matches_reference_env()) {
    expect_no_problems(problems, what)
  } else if (length(problems)) {
    message(sprintf("[informational, environment differs from the stored reference] %s:\n  %s",
                    what, paste(problems, collapse = "\n  ")))
  }
}

# ---- inputs --------------------------------------------------------------------------------

invitro_input <- function() {
  d <- read_fx("inputs", "BIDdetect_data_invitro_delpos.txt.gz")
  # D14: the sample sheet labels noPUS samples treat = input (the legacy script recodes them itself)
  d$treat <- ifelse(d$vector == "noPUS" | d$treat %in% c("in", "input"), "input", "BS")
  d
}

# Deterministic in-cellulo subset, chosen from the stored references (fixed files):
#   1. every site whose WT_mod_Both category differs between the shipped paper table and the
#      Sherlock rerun (paper-vs-rerun drift, e.g. RPL22_chr1_6186768);
#   2. per WT_mod condition: the 4 Modified sites with the largest delta (positive path; the legacy
#      quartile plots need >= 4), and for each other category the 2 sites whose p.adjust_diff is
#      closest to fdr = 0.05 on the log scale (near the decision boundary);
#   3. per PUS7-dependency run: the 2 most and the 2 least significant called sites (the last
#      ones sit at the cutoff), so every run has called sites;
#   4. 12 further oligos drawn with set.seed(20261105) from the rest.
# Per-site results do not depend on which other sites are present: every model is fitted per
# site, the legacy "BH" is a per-site no-op (R-30) and the port runs with p_adjust = "legacy".
# (Only the quartile tables use ntile() across sites; they are not compared.)
incell_subset_chrs <- function(all_chr) {
  rr <- function(f) read_fx("legacy_rerun", "incell", f)
  shipped <- read_fx("shipped", "WT_mod_Both_significance.tsv")
  both <- rr("WT_mod_Both_significance.tsv")
  drift <- dplyr::inner_join(shipped, both, by = c("chr", "pos"), suffix = c(".paper", ".rerun"))
  pick <- drift$chr[drift$category.paper != drift$category.rerun]
  for (cond in c("HepG2", "293T", "Both")) {
    w <- rr(paste0("WT_mod_", cond, "_significance.tsv"))
    mod <- w[w$category == "Modified", ]
    pick <- c(pick, utils::head(mod$chr[order(-mod$delta_delrate, mod$chr)], 4))
    for (cat in setdiff(sort(unique(w$category)), "Modified")) {
      x <- w[w$category == cat, ]
      dist <- abs(log10(pmax(x$p.adjust_diff, 1e-300) / 0.05))
      dist[is.na(dist)] <- Inf
      pick <- c(pick, utils::head(x$chr[order(dist, x$chr)], 2))
    }
  }
  for (f in sort(Sys.glob(fx("legacy_rerun", "incell", "PUS7_dep_*_WT_v_*_significant_summary.tsv")))) {
    s <- read_out(f)
    s <- s[order(s$p.value, s$chr), ]
    pick <- c(pick, utils::head(s$chr, 2), utils::tail(s$chr, 2))
  }
  pick <- unique(pick)
  rest <- sort(setdiff(all_chr, pick))
  set.seed(20261105)
  sort(unique(c(pick, sample(rest, 12))))
}

incell_input <- function() {
  d <- read_fx("inputs", "BIDdetect_data_incell_delpos.txt.gz")
  d$treat <- ifelse(d$treat %in% c("in", "input"), "input", "BS")
  if (!full_run()) d <- d[d$chr %in% incell_subset_chrs(unique(d$chr)), ]
  d
}

# The legacy in-cellulo run is shared by 11b and 11c (it does both analyses in one pass).
legacy_incell <- local({
  cache <- NULL
  function(d) {
    if (is.null(cache)) {
      tmp <- tempfile("legacy_incell_")
      dir.create(tmp)
      input <- file.path(tmp, "input.tsv")
      readr::write_tsv(d, input)
      script <- configure_incell(input, file.path(tmp, "out"), file.path(tmp, "incell_analysis.R"))
      run_legacy(script, wd = tmp)
      cache <<- file.path(tmp, "out", "data_summary")
    }
    cache
  }
})

# ---- 11a: in vitro --------------------------------------------------------------------------

test_that("11a: in vitro treatment test equals modification_analysis.R exactly", {
  skip_if_not(file.exists(fx("inputs", "BIDdetect_data_invitro_delpos.txt.gz")), "fixtures missing")
  skip_without_legacy_pkgs()
  tmp <- tempfile("legacy_invitro_")
  dir.create(tmp)
  run_legacy(normalizePath(fx("legacy_scripts", "modification_analysis.R")),
             c("--input", shQuote(normalizePath(fx("inputs", "BIDdetect_data_invitro_delpos.txt.gz"))),
               "--outdir", shQuote(tmp), "--prefix", "invitro_delpos", "--cores", "1"), wd = tmp)
  legacy_file <- list.files(tmp, "invitro_delpos_modification_significance\\.tsv$", recursive = TRUE, full.names = TRUE)
  expect_length(legacy_file, 1)
  legacy <- read_out(legacy_file)

  r <- suppressWarnings(treatment_test(invitro_input(), p_adjust = "legacy"))
  expect_equal(nrow(r), 770)
  expect_no_problems(diff_invitro(r, legacy), "port vs modification_analysis.R (same environment)")

  check_stored(diff_invitro(r, read_fx("legacy_rerun", "invitro", "invitro_delpos_modification_significance.tsv")),
               "port vs stored Sherlock rerun (in vitro)")
  shipped <- read_fx("shipped", "invitro_modification_significance.tsv")
  m <- dplyr::inner_join(r, shipped, by = c("chr", "pos"), suffix = c("", ".paper"))
  check_stored(c(if (nrow(m) != nrow(shipped)) "rows differ",
                 same_problem(m$category, m$category.paper, "category")),
               "port vs shipped paper table (in vitro categories)")
})

# ---- 11b: in-cellulo WT modification --------------------------------------------------------

wt_mod_runs <- list(
  HepG2 = list(celltype = "HepG2", re = "(1|vector)"),
  `293T` = list(celltype = "293T", re = "(1|vector)"),
  Both = list(celltype = c("HepG2", "293T"), re = "(1|vector) + (1|celltype)"))

test_that("11b: in-cellulo WT_mod treatment tests equal incell_analysis.R", {
  skip_if_not(file.exists(fx("inputs", "BIDdetect_data_incell_delpos.txt.gz")), "fixtures missing")
  skip_without_legacy_pkgs()
  d <- incell_input()
  out <- legacy_incell(d)
  for (cond in names(wt_mod_runs)) {
    run <- wt_mod_runs[[cond]]
    ds <- subset_rows(d, list(celltype = run$celltype, vector = c("WT", "P102", "P4")))
    r <- suppressWarnings(treatment_test(ds, random_effects = run$re, p_adjust = "legacy"))
    legacy <- read_out(file.path(out, paste0("WT_mod_", cond, "_significance.tsv")))
    expect_no_problems(diff_wt_mod(r, legacy, cond), paste("port vs incell_analysis.R (same environment),", cond))

    stored <- read_fx("legacy_rerun", "incell", paste0("WT_mod_", cond, "_significance.tsv")) |>
      dplyr::semi_join(r, by = c("chr", "pos"))
    check_stored(diff_wt_mod(r, stored, cond), paste("port vs stored Sherlock rerun,", cond))
    if (cond == "Both" && full_run()) {
      # Paper vs the stored rerun of the same legacy code: one site changes category because the
      # original fit differed (p 1.4e-5 vs 3.0e-22; package/optimizer drift, not the port).
      shipped <- read_fx("shipped", "WT_mod_Both_significance.tsv")
      mm <- dplyr::inner_join(r, shipped, by = c("chr", "pos"), suffix = c("", ".paper"))
      changed <- mm[mm$category != mm$category.paper, ]
      check_stored(same_problem(sort(paste0(changed$chr, ":", changed$category.paper, ">", changed$category)),
                                "RPL22_chr1_6186768:Inconclusive>Unmodified", "paper drift"),
                   "port vs shipped paper table (WT_mod_Both categories)")
    }
  }
})

# ---- 11c: PUS7 dependency and union ---------------------------------------------------------

pus7_runs <- list(
  HepG2_WT_v_KD = list(celltype = "HepG2", levels = c("P101", "P102"), re = ""),
  HepG2_WT_v_OE = list(celltype = "HepG2", levels = c("P4", "P3"), re = ""),
  `293T_WT_v_KD` = list(celltype = "293T", levels = c("P101", "P102"), re = ""),
  `293T_WT_v_OE` = list(celltype = "293T", levels = c("P4", "P3"), re = ""),
  Both_WT_v_KD = list(celltype = c("HepG2", "293T"), levels = c("P101", "P102"), re = "(1|celltype)"),
  Both_WT_v_OE = list(celltype = c("HepG2", "293T"), levels = c("P4", "P3"), re = "(1|celltype)"))

test_that("11c: PUS7-dependency factor tests and their union equal incell_analysis.R", {
  skip_if_not(file.exists(fx("inputs", "BIDdetect_data_incell_delpos.txt.gz")), "fixtures missing")
  skip_without_legacy_pkgs()
  d <- incell_input()
  out <- legacy_incell(d)
  sites <- dplyr::distinct(d, chr, pos)
  sig <- list()
  for (run_name in names(pus7_runs)) {
    run <- pus7_runs[[run_name]]
    ds <- subset_rows(d, list(celltype = run$celltype, vector = run$levels))
    f <- suppressWarnings(factor_test(ds, factor = "vector", levels = run$levels, random_effects = run$re,
                                      p_adjust = "legacy"))
    sig[[run_name]] <- f$significant
    legacy_file <- file.path(out, paste0("PUS7_dep_", run_name, "_significant_summary.tsv"))
    legacy <- if (file.exists(legacy_file)) read_out(legacy_file) else f$significant[0, ]
    expect_no_problems(diff_factor(f$significant, legacy, run_name, run$levels),
                       paste("port vs incell_analysis.R (same environment),", run_name))

    stored_file <- fx("legacy_rerun", "incell", paste0("PUS7_dep_", run_name, "_significant_summary.tsv"))
    stored <- if (file.exists(stored_file)) read_out(stored_file) else f$significant[0, ]
    check_stored(diff_factor(f$significant, dplyr::semi_join(stored, sites, by = c("chr", "pos")), run_name, run$levels),
                 paste("port vs stored Sherlock rerun,", run_name))
  }
  u <- site_sets(sig, "union")
  legacy_u <- read_out(file.path(out, "PUS7_dep_union_significant_summary.tsv"))
  expect_no_problems(diff_union(u, legacy_u, names(pus7_runs)), "port vs incell_analysis.R (same environment), union")

  stored_u <- read_fx("legacy_rerun", "incell", "PUS7_dep_union_significant_summary.tsv") |>
    dplyr::semi_join(sites, by = c("chr", "pos"))
  check_stored(diff_union(u, stored_u, names(pus7_runs)), "port vs stored Sherlock rerun, union")
  shipped_u <- read_fx("shipped", "PUS7_dep_union_significant_summary.tsv") |> dplyr::semi_join(sites, by = c("chr", "pos"))
  check_stored(same_problem(sort(paste(u$chr, u$pos)), sort(paste(shipped_u$chr, shipped_u$pos)), "union sites"),
               "port vs shipped paper union")
})
