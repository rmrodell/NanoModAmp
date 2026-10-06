# Read funnel (plan §5.2)

`metrics/read_funnel.tsv`, tab-separated, one row per sample and executed step.

Columns: `sample_id, step, file, records, bytes`.

Steps, in this order:

1. `raw` (`.gz` counted correctly; R-08)
2. `trim1_sense`, `trim1_antisense`, `antisense_rc`, `merged` — or `trim1_3prime` instead of these four with `ont_adapter_mode=three_prime_only` (D31)
3. `umi_extracted`
4. `trim2_pool` (MPRA with `orientation_adapters=ont` only)
5. `mapped`, `primary`, `mapq_filtered`
6. `dedup`

Steps a run's parameters skip are **omitted, not written as zero**: `umi_extracted` and `dedup`
when `umi=false`; `trim2_pool` when `orientation_adapters=pool` or `library_type=endogenous`;
the four orientation steps when `ont_adapter_mode=three_prime_only`.

`preprocess/<sample>/<sample>.both_orientations.tsv` records `n_reads_in_both_orientations`
(read IDs found in both the sense and antisense outputs; not produced in `three_prime_only` mode).

Example: [examples/read_funnel.tsv](examples/read_funnel.tsv).
