# Contracts (plan §5)

These files are the interfaces between work packages. They are **frozen after gate G0**; only
the Orchestrator changes them, and every change is logged in [CHANGELOG.md](CHANGELOG.md).

| Contract | File | Machine-readable schema | Example |
| -------- | ---- | ----------------------- | ------- |
| Sample sheet (§5.1) | [samplesheet.md](samplesheet.md) | `assets/schema_input.json` | [examples/samplesheet.csv](examples/samplesheet.csv) |
| Read funnel (§5.2) | [read_funnel.md](read_funnel.md) | — | [examples/read_funnel.tsv](examples/read_funnel.tsv) |
| Counts (§5.3) | [counts.md](counts.md) | — | [examples/sample.counts.tsv](examples/sample.counts.tsv), [examples/counts_merged.tsv](examples/counts_merged.tsv) |
| Analyses config (§5.4) | [analyses.md](analyses.md) | `assets/schema_analyses.json` | [examples/analyses.yaml](examples/analyses.yaml) |
| Site-calling tables (§5.5) | [site_calls.md](site_calls.md) | — | [examples/](examples/) `*_modification_significance.tsv`, `*_significant_summary.tsv`, `*_all_tests.tsv` |
| Parameters (§3) | — | `nextflow_schema.json` | — |
