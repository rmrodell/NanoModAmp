# Stub test data

Tiny placeholder inputs for `-profile test -stub` runs (WP0). They exercise the wiring only:
file and directory FASTQ inputs (`S2_input` is a directory of `*.fastq.gz`), relative paths
resolved against the sample sheet, metadata columns, a treatment and a factor analysis, and
a site set. The synthetic dataset from WP1 replaces them in `conf/test.config`.

`counts/run_a/` and `counts/run_b/` hold two small `counts_merged.tsv` tables (§5.3) for the
`--input_counts` path (D32): merging count tables from runs that were preprocessed separately.

`targets.bed` is standard 0-based (`bed0`, the default). `targets_one_based.bed` is the same sites in the
legacy single-site style (start = end = 1-based site) for `--bed_coordinates one_based_start` tests (R-15).
