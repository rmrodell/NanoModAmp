# Site calling (plan §6.3, §6.4; contracts §5.4, §5.5).
#
# Ported from PUS7regulation2026 @ 4fd285b:
#   Figure3/mpra_sites/modification_analysis.R::diff_equiv_analysis (treatment test, TOST labels)
#   Figure3/mpra_sites/incell_analysis.R::{diff_equiv_analysis, remove_low_sites,
#     perform_anova_test, run_factor_dependency_analysis} (random effects, factor test)
# Model formulas, filters, BH sets and labels are kept exactly as written (rule 0.1, D13, D15–D18,
# D25). Approved changes: seeded future_map with workers = task.cpus (R-03); `na.rm = TRUE` in
# the all-below-SESOI pre-filter (R-17); explicit `treat` values only (D14, no noPUS recoding).
#
# Multiple testing (R-30, PROPOSED): in both legacy scripts `p.adjust(..., "BH")` runs inside a
# mutate() on a tibble still grouped by (chr, pos), so BH sees one p-value at a time and the
# "adjusted" column equals the raw p-value (confirmed in the shipped Figure 3 tables).
# `p_adjust = "legacy"` (default, reproduces the paper) keeps that; `p_adjust = "BH"` applies BH
# over all tested sites as plan §6.3 step 5 / §6.4 step 3 describe. Becca decides at a gate.

# BH over one analysis: legacy = per-site (identity), BH = across all tested sites.
adjust_p <- function(p, p_adjust = c("legacy", "BH")) {
  p_adjust <- match.arg(p_adjust)
  if (p_adjust == "legacy") {
    vapply(p, function(x) stats::p.adjust(x, method = "BH"), numeric(1))
  } else {
    stats::p.adjust(p, method = "BH")
  }
}

#' Validate and normalise a counts table for calling.
#'
#' Requires `chr, pos, treat, rep, delrate, totalReads`; `treat` must be `input` or `BS` (D14).
#' @param data data.frame with the counts_merged.tsv columns (§5.3).
#' @return the same rows as a tibble.
#' @export
check_calling_input <- function(data) {
  required <- c("chr", "pos", "treat", "rep", "delrate", "totalReads")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop("counts table is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  bad <- setdiff(unique(as.character(data$treat)), c("input", "BS"))
  if (length(bad) > 0) {
    stop("treat must be 'input' or 'BS' (D14); found: ", paste(bad, collapse = ", "), call. = FALSE)
  }
  dplyr::as_tibble(data)
}

#' Keep rows whose columns match the allowed values of an analysis `subset` (§5.4).
#' @param data counts table.
#' @param subset named list: column -> allowed values. NULL or empty keeps everything.
#' @export
subset_rows <- function(data, subset = NULL) {
  if (is.null(subset) || length(subset) == 0) return(data)
  for (col in names(subset)) {
    if (!col %in% names(data)) {
      stop("subset column '", col, "' is not in the counts table", call. = FALSE)
    }
    data <- data[as.character(data[[col]]) %in% as.character(unlist(subset[[col]])), , drop = FALSE]
  }
  data
}

# "(1|vector) + (1|celltype)" -> "+ (1|vector) + (1|celltype)", "" -> "" (the legacy scripts pass
# the extra terms with a leading "+"; the analyses YAML stores them without it).
re_suffix <- function(random_effects) {
  re <- trimws(if (is.null(random_effects)) "" else random_effects)
  if (re == "") return("")
  if (!startsWith(re, "+")) re <- paste("+", re)
  re
}

# Run `f` over a list, sequentially or with seeded future_map (R-03). Model fits use no RNG, so
# the result does not depend on the number of workers (testthat 8).
map_sites <- function(x, f, workers = 1) {
  if (workers > 1) {
    old <- future::plan(future::multisession, workers = workers)
    on.exit(future::plan(old), add = TRUE)
  } else {
    old <- future::plan(future::sequential)
    on.exit(future::plan(old), add = TRUE)
  }
  furrr::future_map(x, f, .options = furrr::furrr_options(seed = TRUE))
}

# Equivalence label from the 90% Wald CI of treatBS (modification_analysis.R case_when).
tost_label <- function(lower_ci, upper_ci, eqbound) {
  is_equiv <- (lower_ci > -eqbound) && (upper_ci < eqbound)
  dplyr::case_when(
    is_equiv ~ "Equivalent by TOST",
    lower_ci > eqbound ~ "Different (Positive)",
    upper_ci < -eqbound ~ "Different (Negative)",
    upper_ci > eqbound & lower_ci > -eqbound ~ "Inconclusive (near threshold)",
    lower_ci < -eqbound & upper_ci < eqbound ~ "Unmodified",
    TRUE ~ "Inconclusive (high variance / low power)"
  )
}

# One site: bglmer treatment model, Type III Wald p-value, TOST (§6.3 step 4).
fit_treatment_site <- function(nested_data, model_formula, sesoi) {
  out <- dplyr::tibble(p.value = NA_real_, is_equivalent = NA,
                       equivalence_status = NA_character_, site_specific_eqbound = NA_real_)
  nested_data$treat <- factor(nested_data$treat, levels = c("input", "BS"))
  if (nrow(nested_data) < 4 || length(unique(nested_data$treat)) < 2 ||
      length(unique(nested_data$rep)) < 2) {
    out$equivalence_status <- "Not enough data"
    return(out)
  }
  glm_model <- tryCatch(
    blme::bglmer(model_formula, data = nested_data, family = "binomial",
                 weights = totalReads, fixef.prior = normal),
    error = function(e) NULL)
  if (is.null(glm_model)) {
    out$equivalence_status <- "Model fitting failed"
    return(out)
  }
  p_anova <- tryCatch(car::Anova(glm_model, type = "III"), error = function(e) NULL)
  p_val <- if (!is.null(p_anova) && "treat" %in% rownames(p_anova)) p_anova["treat", "Pr(>Chisq)"] else NA_real_
  baseline_rate <- mean(nested_data$delrate[nested_data$treat == "input"], na.rm = TRUE)
  # D25: when baseline_rate = 0 the bound is infinite; kept as written.
  eqbound <- if (!is.na(baseline_rate) && baseline_rate >= 0 && (baseline_rate + sesoi) < 1) {
    log((baseline_rate + sesoi) / (1 - (baseline_rate + sesoi))) - log(baseline_rate / (1 - baseline_rate))
  } else NA_real_
  is_equiv <- NA
  note <- "Eq. test not run"
  if (!is.na(eqbound)) {
    ci <- tryCatch(stats::confint(glm_model, parm = "treatBS", method = "Wald", level = 0.90),
                   error = function(e) NULL)
    if (!is.null(ci)) {
      is_equiv <- (ci[1, 1] > -eqbound) && (ci[1, 2] < eqbound)
      note <- tost_label(ci[1, 1], ci[1, 2], eqbound)
    }
  }
  dplyr::tibble(p.value = p_val, is_equivalent = is_equiv,
                equivalence_status = note, site_specific_eqbound = eqbound)
}

#' Per-site summary and all-below-SESOI flag (§6.3 step 2).
#' @export
site_summary <- function(data, sesoi = 0.05) {
  data |>
    dplyr::group_by(chr, pos) |>
    dplyr::summarise(
      avg_delrate_BS = mean(delrate[treat == "BS"], na.rm = TRUE),
      avg_delrate_input = mean(delrate[treat == "input"], na.rm = TRUE),
      all_below_thresh = all(delrate < sesoi, na.rm = TRUE),  # R-17
      .groups = "drop") |>
    dplyr::mutate(delta_delrate = avg_delrate_BS - avg_delrate_input)
}

#' Assign Modified / Unmodified / Inconclusive (§6.3 step 6; first match wins).
#' @param x table with all_below_thresh, p.adjust_diff, delta_delrate, is_equivalent.
#' @export
classify_sites <- function(x, sesoi = 0.05, fdr = 0.05) {
  dplyr::mutate(x, category = dplyr::case_when(
    all_below_thresh == TRUE ~ "Unmodified",
    !is.na(p.adjust_diff) & p.adjust_diff < fdr & delta_delrate > sesoi ~ "Modified",
    !is.na(is_equivalent) & is_equivalent == TRUE ~ "Unmodified",
    TRUE ~ "Inconclusive"
  ))
}

#' Treatment analysis: BS vs input per site (§6.3).
#'
#' Model `delrate ~ treat + (1|rep) [+ random_effects]`, binomial bglmer weighted by totalReads
#' with a normal fixed-effect prior; Type III Wald p-value for `treat`; BH over tested sites
#' only (see `p_adjust`); TOST on the 90% Wald CI of treatBS with the site-specific logit bound.
#'
#' @param data counts table (already subset to the analysis).
#' @param sesoi smallest effect size of interest (pre-filter, TOST bound, Δ filter).
#' @param fdr BH threshold for Modified.
#' @param random_effects extra terms, e.g. "(1|vector) + (1|celltype)".
#' @param workers number of parallel workers (task.cpus).
#' @param p_adjust "legacy" (paper: per-site, no correction) or "BH" (across tested sites); R-30.
#' @return tibble in the §5.5 column order, sorted by delta_delrate descending, with
#'   `p.value`, `is_equivalent` and `all_below_thresh` appended.
#' @export
treatment_test <- function(data, sesoi = 0.05, fdr = 0.05, random_effects = "", workers = 1,
                           p_adjust = c("legacy", "BH")) {
  p_adjust <- match.arg(p_adjust)
  data <- check_calling_input(data)
  model_formula <- stats::as.formula(paste("delrate ~ treat + (1|rep)", re_suffix(random_effects)))
  master <- site_summary(data, sesoi)
  low <- dplyr::filter(master, all_below_thresh == TRUE) |> dplyr::select(chr, pos)
  testing <- dplyr::anti_join(data, low, by = c("chr", "pos"))

  nested <- testing |> dplyr::group_by(chr, pos) |> tidyr::nest() |> dplyr::ungroup()
  res <- map_sites(nested$data, function(d) fit_treatment_site(d, model_formula, sesoi), workers)
  stats_tbl <- dplyr::bind_cols(dplyr::select(nested, chr, pos), dplyr::bind_rows(res))
  stats_tbl$p.adjust_diff <- adjust_p(stats_tbl$p.value, p_adjust)

  master |>
    dplyr::left_join(stats_tbl, by = c("chr", "pos")) |>
    classify_sites(sesoi = sesoi, fdr = fdr) |>
    dplyr::select(chr, pos, category, avg_delrate_BS, avg_delrate_input, delta_delrate,
                  p.adjust_diff, equivalence_status, site_specific_eqbound,
                  p.value, is_equivalent, all_below_thresh) |>
    dplyr::arrange(dplyr::desc(delta_delrate))
}

# One site: full vs reduced bglmer, likelihood-ratio test (§6.4 step 2). The model objects must
# be named glm_model_factor / glm_model_nofactor: anova() labels its rows with those names and
# the legacy code keeps the row named "glm_model_factor".
fit_factor_site <- function(nested_data, formula_factor, formula_nofactor, factor_column) {
  if (nrow(nested_data) < 2 || length(unique(nested_data$treat)) < 2 ||
      length(unique(nested_data$rep)) < 2 || length(unique(nested_data[[factor_column]])) < 2) {
    return(dplyr::tibble(term = NA_character_, p.value = NA_real_, note = "Not enough data"))
  }
  glm_model_factor <- tryCatch(
    blme::bglmer(formula_factor, data = nested_data, family = "binomial",
                 weights = totalReads, fixef.prior = normal),
    error = function(e) NULL)
  glm_model_nofactor <- tryCatch(
    blme::bglmer(formula_nofactor, data = nested_data, family = "binomial",
                 weights = totalReads, fixef.prior = normal),
    error = function(e) NULL)
  if (is.null(glm_model_factor) || is.null(glm_model_nofactor)) {
    return(dplyr::tibble(term = NA_character_, p.value = NA_real_, note = "Model fitting failed"))
  }
  p_anova <- tryCatch(stats::anova(glm_model_factor, glm_model_nofactor, test = "ChiSq"),
                      error = function(e) NULL)
  if (is.null(p_anova)) {
    return(dplyr::tibble(term = NA_character_, p.value = NA_real_, note = "ANOVA failed"))
  }
  dplyr::as_tibble(p_anova, rownames = "term") |>
    dplyr::select(term, p.value = `Pr(>Chisq)`) |>
    dplyr::mutate(note = "")
}

#' Factor analysis: does the BS-vs-input effect differ between two levels of a column (§6.4)?
#'
#' Full model `delrate ~ (1|rep) [+RE] + <factor> + L1andBS.ind + L2andBS.ind` vs reduced
#' `delrate ~ (1|rep) [+RE] + <factor> + BS.ind`, likelihood-ratio test, BH over the sites that
#' produced a test. No low-rate pre-filter. Effect: mean delrate per (site, level, treat),
#' `delta_delrate_<level> = BS - input`, `dd_delrate = experimental - baseline`.
#'
#' @param data counts table (already subset to the analysis).
#' @param factor sample-sheet column, e.g. "vector".
#' @param levels c(baseline, experimental).
#' @param direction "positive" (dd > sesoi) or "both" (|dd| > sesoi).
#' @param p_adjust "legacy" (paper) or "BH"; see R-30.
#' @return list(all_tests, significant, baseline_high, experimental_high); the last two only
#'   for direction "both".
#' @export
factor_test <- function(data, factor, levels, random_effects = "", sesoi = 0.05, fdr = 0.05,
                        direction = c("positive", "both"), workers = 1, p_adjust = c("legacy", "BH")) {
  direction <- match.arg(direction)
  p_adjust <- match.arg(p_adjust)
  data <- check_calling_input(data)
  if (!factor %in% names(data)) stop("factor column '", factor, "' is not in the counts table", call. = FALSE)
  if (length(levels) != 2) stop("levels must be [baseline, experimental]", call. = FALSE)
  baseline <- as.character(levels[[1]])
  experimental <- as.character(levels[[2]])
  re <- re_suffix(random_effects)
  formula_factor <- stats::as.formula(paste("delrate ~ (1|rep)", re, "+", factor, "+ L1andBS.ind + L2andBS.ind"))
  formula_nofactor <- stats::as.formula(paste("delrate ~ (1|rep)", re, "+", factor, "+ BS.ind"))

  prepared <- data |>
    dplyr::mutate(!!rlang::sym(factor) := base::factor(as.character(.data[[factor]]), levels = c(baseline, experimental))) |>
    dplyr::mutate(BS.ind = ifelse(treat == "BS", 1, 0),
                  L1andBS.ind = ifelse(treat == "BS" & .data[[factor]] == baseline, 1, 0),
                  L2andBS.ind = ifelse(treat == "BS" & .data[[factor]] == experimental, 1, 0))

  nested <- prepared |> dplyr::group_by(chr, pos) |> tidyr::nest() |> dplyr::ungroup()
  res <- map_sites(nested$data,
                   function(d) fit_factor_site(d, formula_factor, formula_nofactor, factor), workers)
  tests <- dplyr::bind_rows(lapply(seq_along(res), function(i)
    dplyr::bind_cols(nested[i, c("chr", "pos")], res[[i]])))

  kept <- dplyr::filter(tests, term == "glm_model_factor")
  kept$p.value.BH <- adjust_p(kept$p.value, p_adjust)

  deltas <- prepared |>
    dplyr::group_by(chr, pos, .data[[factor]], treat) |>
    dplyr::summarise(mean_delrate = mean(delrate, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = treat, values_from = mean_delrate)
  for (col in c("BS", "input")) if (!col %in% names(deltas)) deltas[[col]] <- NA_real_
  deltas <- deltas |>
    dplyr::mutate(delta_delrate = BS - input) |>
    dplyr::select(chr, pos, dplyr::all_of(factor), delta_delrate) |>
    tidyr::pivot_wider(names_from = dplyr::all_of(factor), values_from = delta_delrate,
                       names_prefix = "delta_delrate_", names_expand = TRUE)
  dcols <- paste0("delta_delrate_", c(baseline, experimental))
  for (col in dcols) if (!col %in% names(deltas)) deltas[[col]] <- NA_real_
  deltas <- deltas |>
    dplyr::select(chr, pos, dplyr::all_of(dcols)) |>
    dplyr::mutate(dd_delrate = .data[[dcols[2]]] - .data[[dcols[1]]])

  cols <- c("chr", "pos", "term", "p.value", "note", "p.value.BH", dcols, "dd_delrate")
  failed <- dplyr::filter(tests, is.na(term) | term != "glm_model_factor")
  failed <- dplyr::filter(failed, !paste(chr, pos) %in% paste(kept$chr, kept$pos))
  all_tests <- dplyr::bind_rows(kept, failed) |>
    dplyr::left_join(deltas, by = c("chr", "pos")) |>
    dplyr::select(dplyr::all_of(cols)) |>
    dplyr::arrange(p.value.BH, chr, pos)

  keep_dir <- if (direction == "positive") function(dd) dd > sesoi else function(dd) abs(dd) > sesoi
  significant <- kept |>
    dplyr::filter(p.value.BH < fdr) |>
    dplyr::left_join(deltas, by = c("chr", "pos")) |>
    dplyr::filter(keep_dir(dd_delrate)) |>
    dplyr::select(dplyr::all_of(cols))

  out <- list(all_tests = all_tests, significant = significant)
  if (direction == "both") {
    out$baseline_high <- dplyr::filter(significant, dd_delrate < 0)
    out$experimental_high <- dplyr::filter(significant, dd_delrate > 0)
  }
  out
}
