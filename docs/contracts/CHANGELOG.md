# Contracts changelog

| Date | Contract | Change | By |
|---|---|---|---|
| 2026-10-06 | all | Initial version from plan §3 and §5 (WP0), including D29 (`orientation_adapters`, `pool_adapter_antisense_5p/3p`, `umi`) and D31 (`ont_adapter_mode`) and the corresponding read-funnel step rules. | Orchestrator |
| 2026-10-06 | counts (§5.3), parameters (§3) | D32: `input_counts` (merge `counts_merged.tsv` tables from earlier runs, then call); new output `counts/merge_sources.tsv`; `--input` and `--input_counts` mutually exclusive. | Orchestrator |
