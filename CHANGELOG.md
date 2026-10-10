# rmrodell/nanomodamp: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased<!-- TODO nf-core: replace with date on release -->]

Initial release of rmrodell/nanomodamp, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- WP0 scaffold: parameters (plan §3) including the opt-in legacy-library modes `orientation_adapters`, `umi` (D29) and `ont_adapter_mode` (D31); sample sheet and analyses schemas; contracts (§5); stubbed preprocessing, counting and calling subworkflows; `test` (stub data) and `test_golden` (placeholder) profiles; CI workflows.
- `--bed_coordinates` defaults to `bed0` (D27, R-15); `one_based_start` reads the paper's legacy BEDs; start = end rows under `bed0` fail validation.
- `--input_counts` (D32): merge `counts_merged.tsv` tables from separately preprocessed runs and call sites; writes `counts/merge_sources.tsv`. `test_golden` endogenous runs as two counting runs plus one merged calling run.
- WP2 preprocessing (plan §6.1): real CAT_FASTQ, linked cutadapt passes with both adapters required (D4), RC + merge with a reads-in-both-orientations count, umi_tools extract, MPRA pool trim, minimap2 (args by library type), `-F 2304` / `-q 30` filters, UMICollapse, and the read funnel (§5.2); the opt-in `pool` / `umi=false` (D29) and `three_prime_only` (D31) modes; fixtures and nf-test specs.

### `Fixed`

### `Dependencies`

- cutadapt 5.2, seqtk 1.4, umi_tools (nf-core `umitools/extract`), minimap2 + samtools (nf-core `minimap2/align`), samtools 1.24, UMICollapse 1.1.0.

### `Deprecated`
