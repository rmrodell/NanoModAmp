# Site sets (plan §6.5), ported from incell_analysis.R (PUS7 union and run_quartile_analysis).

#' Combine per-analysis site tables (§6.5).
#'
#' `union` is a full join on (chr, pos) of the tables in `of` order, every other column suffixed
#' `_<analysis>` (the legacy PUS7 union); `intersection` an inner join; `difference` the sites of
#' the first table that are in none of the others (columns of the first table, suffixed).
#'
#' @param tables named list: analysis name -> site table with chr, pos.
#' @param op "union", "intersection" or "difference".
#' @param of analysis names, in order.
#' @export
site_sets <- function(tables, op = c("union", "intersection", "difference"), of = names(tables)) {
  op <- match.arg(op)
  missing <- setdiff(of, names(tables))
  if (length(missing) > 0) stop("site set refers to unknown analyses: ", paste(missing, collapse = ", "), call. = FALSE)
  renamed <- lapply(of, function(n) {
    t <- dplyr::as_tibble(tables[[n]])
    dplyr::rename_with(t, ~ paste0(.x, "_", n), .cols = -c(chr, pos))
  })
  names(renamed) <- of
  if (op == "union") {
    out <- dplyr::tibble(chr = character(), pos = integer())
    for (t in renamed) out <- dplyr::full_join(out, t, by = c("chr", "pos"))
  } else if (op == "intersection") {
    out <- Reduce(function(a, b) dplyr::inner_join(a, b, by = c("chr", "pos")), renamed)
  } else {
    out <- renamed[[1]]
    for (t in renamed[-1]) out <- dplyr::anti_join(out, t, by = c("chr", "pos"))
  }
  out
}

#' Bin a site set into quartiles of another analysis's column (run_quartile_analysis).
#'
#' Sites without a value in `ref` are dropped; NULL when fewer than 4 sites remain.
#' @param sites site table with chr, pos.
#' @param ref table with chr, pos and `column` (e.g. a treatment table's delta_delrate).
#' @export
site_quartiles <- function(sites, ref, column) {
  if (!column %in% names(ref)) stop("quartiles_by column '", column, "' not in the reference analysis", call. = FALSE)
  q <- sites |>
    dplyr::left_join(dplyr::select(ref, chr, pos, dplyr::all_of(column)), by = c("chr", "pos")) |>
    dplyr::filter(!is.na(.data[[column]])) |>
    dplyr::mutate(quartile = dplyr::ntile(.data[[column]], 4))
  if (nrow(q) < 4) return(NULL)
  q
}
