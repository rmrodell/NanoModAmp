# Analysis runners used by bin/nma_call.R and bin/nma_sitesets.R: read the analyses YAML
# (§5.4), run one analysis, and write the §4 output layout:
#   calling/<name>/{<name>_log.txt, <name>_resolved_config.yaml, data_summary/, data_raw/, plots/}
#   site_sets/<set>/{data_summary/, data_raw/, plots/}

ANALYSIS_KEYS <- c("name", "type", "subset", "random_effects", "sesoi", "fdr", "plot_all_sites",
                   "colors", "factor", "levels", "direction", "p_adjust")

#' Read and validate an analyses YAML (§5.4).
#' @return list(analyses = named list of resolved analyses, site_sets = list).
#' @export
read_analyses <- function(path, plot_all_sites = FALSE, p_adjust = "BH") {
  cfg <- yaml::read_yaml(path)
  if (is.null(cfg$analyses) || length(cfg$analyses) == 0) stop("analyses YAML has no 'analyses'", call. = FALSE)
  out <- list()
  for (a in cfg$analyses) {
    bad <- setdiff(names(a), ANALYSIS_KEYS)
    if (length(bad) > 0) stop("analysis '", a$name, "': unknown keys ", paste(bad, collapse = ", "), call. = FALSE)
    if (is.null(a$name) || is.null(a$type)) stop("every analysis needs 'name' and 'type'", call. = FALSE)
    if (!a$type %in% c("treatment", "factor")) stop("analysis '", a$name, "': type must be treatment or factor", call. = FALSE)
    if (a$name %in% names(out)) stop("duplicate analysis name '", a$name, "'", call. = FALSE)
    if (a$type == "factor" && (is.null(a$factor) || length(a$levels) != 2)) {
      stop("factor analysis '", a$name, "' needs 'factor' and two 'levels'", call. = FALSE)
    }
    a$random_effects <- if (is.null(a$random_effects)) "" else a$random_effects
    a$sesoi <- if (is.null(a$sesoi)) 0.05 else a$sesoi
    a$fdr <- if (is.null(a$fdr)) 0.05 else a$fdr
    a$p_adjust <- if (is.null(a$p_adjust)) p_adjust else a$p_adjust  # R-30
    if (!a$p_adjust %in% c("BH", "legacy")) stop("analysis '", a$name, "': p_adjust must be BH or legacy", call. = FALSE)
    a$plot_all_sites <- if (is.null(a$plot_all_sites)) plot_all_sites else a$plot_all_sites
    if (a$type == "factor") a$direction <- if (is.null(a$direction)) "positive" else a$direction
    out[[a$name]] <- a
  }
  for (s in cfg$site_sets) {
    if (is.null(s$name) || is.null(s$op) || length(s$of) < 2) {
      stop("every site set needs 'name', 'op' and at least two analyses in 'of'", call. = FALSE)
    }
    if (!s$op %in% c("union", "intersection", "difference")) {
      stop("site set '", s$name, "': op must be union, intersection or difference", call. = FALSE)
    }
    unknown <- setdiff(s$of, names(out))
    if (length(unknown) > 0) stop("site set '", s$name, "' refers to unknown analyses: ", paste(unknown, collapse = ", "), call. = FALSE)
    # Same rule as the pipeline's startup validation: quartiles_by must name a defined analysis.
    qb <- s$quartiles_by
    if (!is.null(qb)) {
      if (is.null(qb$analysis) || is.null(qb$column)) {
        stop("site set '", s$name, "': quartiles_by needs 'analysis' and 'column'", call. = FALSE)
      }
      if (!qb$analysis %in% names(out)) {
        stop("site set '", s$name, "': quartiles_by.analysis '", qb$analysis, "' is not a defined analysis", call. = FALSE)
      }
    }
  }
  list(analyses = out, site_sets = if (is.null(cfg$site_sets)) list() else cfg$site_sets)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  readr::write_tsv(x, path, na = "NA", progress = FALSE)
}

#' Read counts_merged.tsv (§5.3).
#' @export
read_counts <- function(path) {
  readr::read_tsv(path, show_col_types = FALSE, progress = FALSE, guess_max = 1e6)
}

#' Run one analysis and write its outputs under `outdir/<name>/`.
#' @param counts counts table (all samples). @param a one resolved analysis from read_analyses().
#' @return invisible list with the result tables.
#' @export
run_analysis <- function(counts, a, outdir, workers = 1) {
  name <- a$name
  run_dir <- file.path(outdir, name)
  sum_dir <- file.path(run_dir, "data_summary")
  raw_dir <- file.path(run_dir, "data_raw")
  plot_dir <- file.path(run_dir, "plots")
  for (d in c(sum_dir, raw_dir, plot_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
  log_file <- file.path(run_dir, paste0(name, "_log.txt"))
  log <- function(...) cat(..., "\n", file = log_file, append = TRUE, sep = "")
  cat("", file = log_file)
  yaml::write_yaml(a, file.path(run_dir, paste0(name, "_resolved_config.yaml")))

  data <- check_calling_input(subset_rows(counts, a$subset))
  log("analysis: ", name, " (", a$type, ")")
  log("rows: ", nrow(data), "; sites: ", dplyr::n_distinct(data$chr, data$pos), "; samples: ",
      if ("sample_id" %in% names(data)) dplyr::n_distinct(data$sample_id) else NA)
  log("random effects: ", if (a$random_effects == "") "(none)" else a$random_effects,
      "; sesoi: ", a$sesoi, "; fdr: ", a$fdr, "; p_adjust: ", a$p_adjust, "; workers: ", workers)
  if (nrow(data) == 0) stop("analysis '", name, "': subset selects no rows", call. = FALSE)
  write_tsv(data, file.path(raw_dir, paste0(name, "_standardized_input.tsv")))

  if (a$type == "treatment") {
    res <- treatment_test(data, sesoi = a$sesoi, fdr = a$fdr, random_effects = a$random_effects, workers = workers,
                          p_adjust = a$p_adjust)
    log("pre-filtered (all delrate < sesoi): ", sum(res$all_below_thresh), "; tested: ", sum(!is.na(res$p.value)),
        "; not enough data: ", sum(res$equivalence_status %in% "Not enough data"),
        "; model fitting failed: ", sum(res$equivalence_status %in% "Model fitting failed"))
    tab <- table(factor(res$category, levels = c("Modified", "Unmodified", "Inconclusive")))
    log("categories: ", paste(names(tab), tab, sep = "=", collapse = ", "))
    write_tsv(res, file.path(sum_dir, paste0(name, "_modification_significance.tsv")))
    for (cat_name in c("Modified", "Unmodified")) {
      sites <- dplyr::filter(res, category == cat_name)
      tag <- tolower(cat_name)
      write_tsv(sites, file.path(sum_dir, paste0(name, "_", tag, "_sites_summary.tsv")))
      write_tsv(dplyr::semi_join(data, sites, by = c("chr", "pos")),
                file.path(raw_dir, paste0(name, "_", tag, "_sites_raw_data.tsv")))
    }
    modified <- dplyr::filter(res, category == "Modified")
    if (nrow(modified) >= 4) {
      q <- dplyr::mutate(modified, quartile = dplyr::ntile(delta_delrate, 4))
      for (k in sort(unique(q$quartile))) {
        qs <- dplyr::filter(q, quartile == k)
        write_tsv(qs, file.path(sum_dir, paste0(name, "_modified_sites_Q", k, "_summary.tsv")))
        write_tsv(dplyr::semi_join(data, qs, by = c("chr", "pos")),
                  file.path(raw_dir, paste0(name, "_modified_sites_Q", k, "_raw.tsv")))
      }
    } else {
      log("quartile analysis skipped: ", nrow(modified), " Modified sites (< 4)")
    }
    plot_treatment(res, data, plot_dir, name, a$colors, a$plot_all_sites)
    out <- list(table = res)
  } else {
    res <- factor_test(data, factor = a$factor, levels = unlist(a$levels), random_effects = a$random_effects,
                       sesoi = a$sesoi, fdr = a$fdr, direction = a$direction, workers = workers,
                       p_adjust = a$p_adjust)
    log("tested: ", sum(res$all_tests$term %in% "glm_model_factor"),
        "; not enough data: ", sum(res$all_tests$note %in% "Not enough data"),
        "; model fitting failed: ", sum(res$all_tests$note %in% "Model fitting failed"),
        "; significant: ", nrow(res$significant))
    write_tsv(res$all_tests, file.path(sum_dir, paste0(name, "_all_tests.tsv")))
    write_tsv(res$significant, file.path(sum_dir, paste0(name, "_significant_summary.tsv")))
    write_tsv(dplyr::semi_join(data, res$significant, by = c("chr", "pos")),
              file.path(raw_dir, paste0(name, "_significant_raw_data.tsv")))
    if (a$direction == "both") {
      lv <- as.character(unlist(a$levels))
      write_tsv(res$baseline_high, file.path(sum_dir, paste0(name, "_", lv[1], "_high_summary.tsv")))
      write_tsv(res$experimental_high, file.path(sum_dir, paste0(name, "_", lv[2], "_high_summary.tsv")))
    }
    plot_factor(res$significant, data, a$factor, unlist(a$levels), plot_dir, name, a$colors, a$plot_all_sites)
    out <- list(table = res$significant, all_tests = res$all_tests)
  }
  log("done")
  invisible(out)
}

# The table a site set uses from an analysis: Modified sites of a treatment analysis, significant
# sites of a factor analysis.
analysis_sites_table <- function(calling_dir, a) {
  sum_dir <- file.path(calling_dir, a$name, "data_summary")
  if (a$type == "treatment") {
    t <- readr::read_tsv(file.path(sum_dir, paste0(a$name, "_modification_significance.tsv")), show_col_types = FALSE)
    list(sites = dplyr::filter(t, category == "Modified"), full = t)
  } else {
    t <- readr::read_tsv(file.path(sum_dir, paste0(a$name, "_significant_summary.tsv")), show_col_types = FALSE)
    list(sites = t, full = t)
  }
}

#' Build every site set of an analyses config from finished analysis directories (§6.5).
#' @param calling_dir directory holding calling/<name>/ results.
#' @param counts optional counts table for the raw-data tables.
#' @export
run_site_sets <- function(cfg, calling_dir, outdir, counts = NULL) {
  for (s in cfg$site_sets) {
    set_dir <- file.path(outdir, s$name)
    tables <- lapply(cfg$analyses[s$of], function(a) analysis_sites_table(calling_dir, a)$sites)
    names(tables) <- s$of
    res <- site_sets(tables, op = s$op, of = s$of)
    write_tsv(res, file.path(set_dir, "data_summary", paste0(s$name, "_summary.tsv")))
    if (!is.null(counts)) {
      write_tsv(dplyr::semi_join(counts, res, by = c("chr", "pos")),
                file.path(set_dir, "data_raw", paste0(s$name, "_raw_data.tsv")))
    }
    qb <- s$quartiles_by
    if (!is.null(qb)) {
      ref <- analysis_sites_table(calling_dir, cfg$analyses[[qb$analysis]])$full
      q <- site_quartiles(res, ref, qb$column)
      qname <- paste0(s$name, "_by_", qb$analysis)
      if (is.null(q)) {
        cat("site set ", s$name, ": fewer than 4 sites with ", qb$column, "; quartiles skipped\n", sep = "")
      } else {
        write_tsv(q, file.path(set_dir, "data_summary", paste0(qname, "_quartiles.tsv")))
        plot_quartiles(q, qb$column, file.path(set_dir, "plots"), qname)
      }
    }
  }
  invisible(TRUE)
}
