# nanomodamp-r: R image for the nanomodamp R package, the bin/nma_*.R CLIs and the vendored legacy
# site-calling scripts (plan §2, D23, R-24). Published to ghcr.io/rmrodell/nanomodamp-r by
# .github/workflows/container.yml.
#
# Every R package is pinned to the version of the Sherlock reference environment (R 4.3.2 +
# Bioconductor 3.18) that the WP3/WP4 tests and the test-11 legacy references were produced with:
# r-packages.tsv (written by make_lockfile.R). The build fails if any version differs.
#
# Build from the repo root:  docker build -f containers/nanomodamp-r/Dockerfile -t nanomodamp-r:dev .
FROM rocker/r-ver:4.3.2@sha256:4e32addfc4da3e660f6e0d05ce5e43d3eceb9db58a60b9a142e0dde9a654ead1

LABEL org.opencontainers.image.source="https://github.com/rmrodell/NanoModAmp" \
      org.opencontainers.image.description="R 4.3.2 / Bioconductor 3.18 environment for the NanoModAmp pipeline (nanomodamp R package)" \
      org.opencontainers.image.licenses="MIT"

# System libraries for the pinned packages (curl/RCurl, xml2, Rhtslib/Rsamtools, png/jpeg,
# nloptr -> cmake, argparse -> python3) and procps (Nextflow uses ps for task metrics).
RUN apt-get update -qq \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      cmake libbz2-dev libcurl4-openssl-dev libdeflate-dev libfontconfig1-dev libfreetype6-dev \
      libfribidi-dev libharfbuzz-dev libjpeg-dev liblzma-dev libpng-dev libssl-dev libtiff-dev \
      libxml2-dev make procps python3 zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*

COPY containers/nanomodamp-r/r-packages.tsv containers/nanomodamp-r/install_lockfile.R /opt/nanomodamp-r/
RUN MAKEFLAGS="-j$(nproc)" Rscript /opt/nanomodamp-r/install_lockfile.R /opt/nanomodamp-r/r-packages.tsv

COPY R/nanomodamp /tmp/nanomodamp
RUN R CMD INSTALL /tmp/nanomodamp && rm -rf /tmp/nanomodamp \
 && Rscript -e 'library(nanomodamp); cat("nanomodamp", as.character(packageVersion("nanomodamp")), "\n")'
