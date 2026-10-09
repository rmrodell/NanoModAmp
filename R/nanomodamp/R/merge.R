# Merging count tables by column name, never by position (plan §5.3, R-06).
# All values are read and written as text, so numbers keep their exact formatting.

read_text_table <- function(file, sep = "\t") {
  fread(file, sep = sep, colClasses = "character", na.strings = c("NA", ""), quote = "")
}

#' Merge per-sample counts into counts_merged.tsv (§5.3).
#'
#' @param count_files Per-sample `<sample_id>.counts.tsv` files (the sample_id is the file name).
#' @param samplesheet Sample sheet CSV (§5.1); its columns other than `fastq` are the metadata.
#' @return data.table: `sample_id`, the metadata columns in sample-sheet order, then the §5.3
#'   count columns. Samples are in sample-sheet order; rows within a sample keep file order.
merge_counts <- function(count_files, samplesheet) {
  ss <- fread(samplesheet, sep = ",", colClasses = "character", na.strings = "")
  if (!"sample_id" %in% names(ss)) stop("sample sheet has no sample_id column", call. = FALSE)
  meta_cols <- setdiff(names(ss), c("sample_id", "fastq"))
  ids <- sub("\\.counts\\.tsv$", "", basename(count_files))
  unknown <- setdiff(ids, ss$sample_id)
  if (length(unknown)) stop("count files for samples not in the sample sheet: ",
                            paste(unknown, collapse = ", "), call. = FALSE)
  dup <- ids[duplicated(ids)]
  if (length(dup)) stop("more than one count file for: ", paste(unique(dup), collapse = ", "), call. = FALSE)
  missing <- setdiff(ss$sample_id, ids)
  if (length(missing)) warning("no count file for sample(s): ", paste(missing, collapse = ", "), call. = FALSE)

  tabs <- lapply(seq_along(count_files), function(i) {
    t <- read_text_table(count_files[i])
    absent <- setdiff(COUNT_COLUMNS, names(t))
    if (length(absent)) stop(basename(count_files[i]), " lacks column(s): ", paste(absent, collapse = ", "), call. = FALSE)
    t <- t[, ..COUNT_COLUMNS]
    t[, sample_id := ids[i]]
    t
  })
  out <- rbindlist(tabs[order(match(ids, ss$sample_id))], use.names = TRUE)
  if (nrow(out) == 0) out <- data.table(sample_id = character())
  if (length(meta_cols)) out[ss, on = "sample_id", (meta_cols) := mget(paste0("i.", meta_cols))]
  cols <- c("sample_id", meta_cols, COUNT_COLUMNS)
  for (c in setdiff(cols, names(out))) out[, (c) := character()]
  out[, ..cols]
}

#' Merge counts_merged.tsv tables from separate runs (D32, `--input_counts`).
#'
#' @param tables Paths of counts_merged.tsv tables.
#' @param sources Labels recorded in merge_sources.tsv (default: the paths).
#' @return list: `merged` (§5.3 layout; metadata = union in first-seen order, missing -> NA) and
#'   `sources` (sample_id, source_table).
merge_count_tables <- function(tables, sources = tables) {
  if (length(tables) == 0) stop("no tables to merge", call. = FALSE)
  if (length(sources) != length(tables)) stop("sources must match tables", call. = FALSE)
  tabs <- lapply(tables, read_text_table)
  meta <- character()
  for (i in seq_along(tabs)) {
    n <- names(tabs[[i]])
    if (!"sample_id" %in% n) stop(sources[i], ": no sample_id column", call. = FALSE)
    if (anyDuplicated(n)) stop(sources[i], ": duplicated column names", call. = FALSE)
    absent <- setdiff(COUNT_COLUMNS, n)
    if (length(absent)) {
      stop(sources[i], ": count/site columns differ from the §5.3 layout (missing: ",
           paste(absent, collapse = ", "), ")", call. = FALSE)
    }
    meta <- c(meta, setdiff(n, c("sample_id", COUNT_COLUMNS, meta)))
  }
  seen <- character()
  for (i in seq_along(tabs)) {
    ids <- unique(tabs[[i]]$sample_id)
    dup <- intersect(ids, seen)
    if (length(dup)) stop("sample_id in more than one table: ", paste(dup, collapse = ", "),
                          " (", sources[i], ")", call. = FALSE)
    seen <- c(seen, ids)
    absent <- setdiff(meta, names(tabs[[i]]))
    if (length(absent)) {
      warning(sources[i], " lacks metadata column(s) ", paste(absent, collapse = ", "),
              "; filled with NA", call. = FALSE)
      for (c in absent) tabs[[i]][, (c) := NA_character_]
    }
  }
  cols <- c("sample_id", meta, COUNT_COLUMNS)
  merged <- rbindlist(lapply(tabs, function(t) t[, ..cols]), use.names = TRUE)
  src <- rbindlist(lapply(seq_along(tabs), function(i)
    data.table(sample_id = unique(tabs[[i]]$sample_id), source_table = sources[i])))
  list(merged = merged, sources = src)
}
