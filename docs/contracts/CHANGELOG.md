# Contracts changelog

| Date | Contract | Change | By |
|---|---|---|---|
| 2026-10-06 | all | Initial version from plan §3 and §5 (WP0), including D29 (`orientation_adapters`, `pool_adapter_antisense_5p/3p`, `umi`) and D31 (`ont_adapter_mode`) and the corresponding read-funnel step rules. | Orchestrator |
| 2026-10-06 | counts (§5.3), parameters (§3) | D32: `input_counts` (merge `counts_merged.tsv` tables from earlier runs, then call); new output `counts/merge_sources.tsv`; `--input` and `--input_counts` mutually exclusive. | Orchestrator |
| 2026-10-06 | parameters (§3) | G1-e (D27, R-15): `bed_coordinates` default `bed0`; start = end rows under `bed0` are a validation error. | Orchestrator |
| 2026-10-09 | analyses (§5.4) | Example treatment analysis changed from in vitro to in cellulo (`incellulo`, HepG2 + 293T) by Becca. Example made valid: added the referenced `PUS7_dep_Both_WT_v_OE` factor analysis (levels [P4 = WT, P3 = OE], paper order) and pointed `quartiles_by` at `incellulo` (it named an undefined `WT_mod_Both`). Startup validation now also checks `quartiles_by.analysis`. | Becca / Orchestrator |
