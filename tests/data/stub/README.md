# Stub test data

Tiny placeholder inputs for `-profile test -stub` runs (WP0). They exercise the wiring only:
file and directory FASTQ inputs (`S2_input` is a directory of `*.fastq.gz`), relative paths
resolved against the sample sheet, metadata columns, a treatment and a factor analysis, and
a site set. The synthetic dataset from WP1 replaces them in `conf/test.config`.

`counts/run_a/` and `counts/run_b/` hold two small `counts_merged.tsv` tables (§5.3) for the
`--input_counts` path (D32): merging count tables from runs that were preprocessed separately.
