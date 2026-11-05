# nanomodamp-r container

R image for the pipeline's R steps (COUNT, MERGE_COUNTS, MERGE_COUNT_TABLES, BACKGROUND_QC,
CALL_SITES, SITE_SETS), the `bin/nma_*.R` CLIs, the CI R tests, and the vendored legacy
site-calling scripts used by the legacy-equivalence test.

**Image:** `ghcr.io/rmrodell/nanomodamp-r` (GitHub Container Registry; decision 2026-11-05).
Tags: `<nanomodamp version>` (e.g. `0.1.0`, used by the modules and CI), `sha-<commit>`, `dev`.

## What is pinned

- **Base:** `rocker/r-ver:4.3.2` by digest
  (`sha256:4e32addfc4da3e660f6e0d05ce5e43d3eceb9db58a60b9a142e0dde9a654ead1`), so R itself and
  its recommended packages (Matrix 1.6-1.1, MASS 7.3-60, nlme 3.1-163, …) are R 4.3.2's.
- **Every other R package** at the exact version of the Sherlock reference environment
  (R 4.3.2 module + `~/R/x86_64-pc-linux-gnu-library/4.3`, Bioconductor 3.18) in which the WP3/WP4
  tests and the test-11 legacy references were produced: `r-packages.tsv` (137 rows: 109 CRAN,
  19 Bioconductor, 9 shipped with R), e.g. lme4 1.1-35.5, blme 1.0-6, car 3.1-3, data.table
  1.16.4, Rsamtools 2.18.0.
- Packages are installed from exact source tarballs (CRAN current or Archive; Bioconductor 3.18
  current or Archive), dependencies first, with no dependency resolution. The build **fails** if
  any package (including the shipped recommended ones) ends up at a different version.

Why it matters: the mixed-model fits (`blme::bglmer` on `lme4`) converge to slightly different
estimates on other package versions — on newer CRAN versions p-values move by 1e-4 to 1e-2
(relative) and a borderline site can change category. Results are reproducible only in this image.

## Files

| File | Purpose |
|---|---|
| `Dockerfile` | Image definition. Build from the repo root: `docker build -f containers/nanomodamp-r/Dockerfile -t nanomodamp-r:dev .` |
| `r-packages.tsv` | Lockfile: `package`, `version`, `repository` (`CRAN`, `Bioconductor`, or `R` = shipped with R, checked only), in install order |
| `make_lockfile.R` | Regenerates `r-packages.tsv` from an R installation |
| `install_lockfile.R` | Installs and verifies the lockfile inside the image |

## Rebuilding the lockfile

Only when the reference environment changes on purpose (e.g. a new R version for a release). On
Sherlock:

```bash
unset R_HOME; module load R/4.3.2
/share/software/user/open/R/4.3.2/bin/Rscript containers/nanomodamp-r/make_lockfile.R .
```

Roots: the `Imports` of `R/nanomodamp/DESCRIPTION`, plus testthat, optparse, argparse, furrr,
readr, forcats, purrr, ggtext, future.callr and GenomicAlignments (CLIs and the vendored legacy
scripts). Changing versions changes results; record it in the Change Register and regenerate the
test-11 references in the new image.

## Publishing

`.github/workflows/container.yml` builds the image on pushes to `dev`/`main` that touch
`containers/nanomodamp-r/` or `R/nanomodamp/` (not its tests), and on manual dispatch, and pushes
it to GHCR with the workflow's own `GITHUB_TOKEN`. Pull requests build without pushing.

**Make the package public** after the first publish, so users and Nextflow can pull it without
logging in: GitHub → your profile → **Packages** → **nanomodamp-r** → **Package settings** →
**Change visibility** → **Public**. CI logs in to GHCR, so it works either way.
