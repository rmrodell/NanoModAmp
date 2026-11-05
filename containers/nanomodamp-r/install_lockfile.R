#!/usr/bin/env Rscript
# Install every package of r-packages.tsv at its exact version, in file order (dependencies
# first), from source tarballs: CRAN (current or Archive) and Bioconductor 3.18 (current or
# Archive; data/annotation repo for data packages). Rows with repository "R" ship with R and are
# only checked. Fails if any package ends up at a different version.
#
# Usage: Rscript install_lockfile.R r-packages.tsv

lock <- read.delim(commandArgs(trailingOnly = TRUE)[1], comment.char = "#",
                   colClasses = "character")
cran <- "https://cloud.r-project.org/src/contrib"
bioc <- c("https://bioconductor.org/packages/3.18/bioc/src/contrib",
          "https://bioconductor.org/packages/3.18/data/annotation/src/contrib",
          "https://bioconductor.org/packages/3.18/data/experiment/src/contrib")

candidates <- function(pkg, ver, repo) {
  tb <- sprintf("%s_%s.tar.gz", pkg, ver)
  bases <- if (repo == "CRAN") cran else bioc
  c(file.path(bases, tb), file.path(bases, "Archive", pkg, tb))
}

installed_version <- function(pkg) {
  v <- tryCatch(as.character(packageVersion(pkg, lib.loc = .libPaths())), error = function(e) NA)
  v
}

for (i in seq_len(nrow(lock))) {
  pkg <- lock$package[i]; ver <- lock$version[i]; repo <- lock$repository[i]
  if (repo == "R") next
  if (identical(installed_version(pkg), as.character(package_version(ver)))) next
  dest <- file.path(tempdir(), sprintf("%s_%s.tar.gz", pkg, ver))
  ok <- FALSE
  for (u in candidates(pkg, ver, repo)) {
    if (tryCatch(utils::download.file(u, dest, quiet = TRUE, mode = "wb") == 0, error = function(e) FALSE,
                 warning = function(w) FALSE)) { ok <- TRUE; break }
  }
  if (!ok) stop(sprintf("no source tarball found for %s %s (%s)", pkg, ver, repo))
  message(sprintf("[%d/%d] %s %s", i, nrow(lock), pkg, ver))
  utils::install.packages(dest, repos = NULL, type = "source", dependencies = FALSE,
                          Ncpus = max(1L, parallel::detectCores()))
  unlink(dest)
}

got <- vapply(lock$package, installed_version, character(1))
want <- as.character(package_version(lock$version))
bad <- which(is.na(got) | got != want)
if (length(bad)) {
  stop("version mismatch:\n", paste(sprintf("  %s: want %s, got %s", lock$package[bad], want[bad],
                                            got[bad]), collapse = "\n"))
}
message(sprintf("all %d lockfile packages at the pinned versions", nrow(lock)))
