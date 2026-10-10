# Small simulated counts tables for the site-calling tests (seeded, deterministic).
#
# sim_site(): one site, `reps` batch-paired replicates per group, deletion counts drawn as
# Binomial(n, p) with a replicate-level jitter so (1|rep) has something to absorb.
sim_site <- function(chr, pos, p, n = 500, reps = 1:3, extra = list(), seed = 1) {
  set.seed(seed)
  rows <- list()
  for (g in names(p)) {
    parts <- strsplit(g, ":", fixed = TRUE)[[1]]  # "<treat>" or "<treat>:<level>"
    for (r in reps) {
      pr <- if (p[[g]] == 0) 0 else min(max(p[[g]] + stats::rnorm(1, 0, 0.005), 0), 1)
      del <- stats::rbinom(1, n, pr)
      row <- data.frame(chr = chr, pos = pos, treat = parts[1], rep = r, totalReads = n,
                        Deletion.count = del, delrate = del / n, stringsAsFactors = FALSE)
      if (length(parts) > 1) row$level <- parts[2]
      for (k in names(extra)) row[[k]] <- extra[[k]]
      rows[[length(rows) + 1]] <- row
    }
  }
  do.call(rbind, rows)
}

bind_sites <- function(...) {
  x <- list(...)
  cols <- unique(unlist(lapply(x, names)))
  x <- lapply(x, function(d) { for (c in setdiff(cols, names(d))) d[[c]] <- NA; d[cols] })
  dplyr::as_tibble(do.call(rbind, x))
}

# A treatment dataset with one clear-cut site of each kind.
sim_treatment <- function() {
  bind_sites(
    sim_site("POS", 10, c(input = 0.01, BS = 0.60), n = 800, seed = 1),            # Modified
    sim_site("NULL", 10, c(input = 0.20, BS = 0.20), n = 4000, seed = 2),          # Unmodified (TOST)
    sim_site("LOW", 10, c(input = 0.01, BS = 0.02), n = 500, seed = 3),            # all < SESOI
    sim_site("NOISY", 10, c(input = 0.15, BS = 0.18), n = 22, reps = 1:2, seed = 4) # Inconclusive
  )
}

# A factor dataset: BS effect bigger in level B (positive interaction), equal effects, and one
# site where the effect is bigger in the baseline level A.
sim_factor <- function() {
  bind_sites(
    sim_site("INTER", 10, c("input:A" = 0.02, "BS:A" = 0.05, "input:B" = 0.02, "BS:B" = 0.50), n = 800, seed = 11),
    sim_site("EQUAL", 10, c("input:A" = 0.02, "BS:A" = 0.40, "input:B" = 0.02, "BS:B" = 0.40), n = 800, seed = 12),
    sim_site("REVERSE", 10, c("input:A" = 0.02, "BS:A" = 0.50, "input:B" = 0.02, "BS:B" = 0.05), n = 800, seed = 13)
  )
}
