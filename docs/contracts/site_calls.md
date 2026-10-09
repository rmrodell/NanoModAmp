# Site-calling result tables (plan §5.5)

All tab-separated, under `calling/<name>/data_summary/`.

**Treatment** (`type: treatment`): `<name>_modification_significance.tsv`. Columns:
`chr, pos, category, avg_delrate_BS, avg_delrate_input, delta_delrate, p.adjust_diff, equivalence_status, site_specific_eqbound`,
sorted by `delta_delrate` descending. May append `p.value`, `is_equivalent`, `model_status`.
`category` ∈ {`Modified`, `Unmodified`, `Inconclusive`}; `equivalence_status` uses the labels
of `modification_analysis.R` (`Equivalent by TOST`, `Different (Positive)`, `Different (Negative)`,
`Inconclusive (near threshold)`, `Unmodified`, `Inconclusive (high variance / low power)`, `Eq. test not run`).

**Factor** (`type: factor`): `<name>_significant_summary.tsv` with
`chr, pos, term, p.value, note, p.value.BH, delta_delrate_<baseline>, delta_delrate_<experimental>, dd_delrate`,
and `<name>_all_tests.tsv` (same columns, all sites including non-significant and failed ones with
their `note`). With `direction: both`, also `<name>_<baseline>_high` (dd < 0) and
`<name>_<experimental>_high` (dd > 0).

Examples: [examples/invitro_modification_significance.tsv](examples/invitro_modification_significance.tsv),
[examples/PUS7_dep_Both_WT_v_KD_significant_summary.tsv](examples/PUS7_dep_Both_WT_v_KD_significant_summary.tsv),
[examples/PUS7_dep_Both_WT_v_KD_all_tests.tsv](examples/PUS7_dep_Both_WT_v_KD_all_tests.tsv).
Values in the examples are illustrative, not results.
