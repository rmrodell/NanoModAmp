# nanomodamp-r: R/Bioconductor image for the nanomodamp R package and bin/nma_*.R CLIs (plan §2, D23).
# Build from the repo root:  docker build -f containers/nanomodamp-r/Dockerfile -t nanomodamp-r:dev .
# Versions match the stack the WP3 tests ran on (R 4.3.2, Bioconductor 3.18). WP4 adds its
# site-calling packages (blme, car, furrr, ...) here. Not yet built or pushed (see PROGRESS.md).
FROM bioconductor/bioconductor_docker:RELEASE_3_18

RUN R -q -e 'BiocManager::install(c("Rsamtools", "GenomicRanges", "IRanges", "S4Vectors", \
        "GenomeInfoDb", "Biostrings"), ask = FALSE, update = FALSE)' \
 && R -q -e 'install.packages(c("data.table", "optparse", "future", "future.apply", "ggplot2", \
        "testthat"), repos = "https://cloud.r-project.org")'

COPY R/nanomodamp /tmp/nanomodamp
RUN R CMD INSTALL /tmp/nanomodamp && rm -rf /tmp/nanomodamp
