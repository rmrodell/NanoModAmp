# Synthetic data simulator (WP1)

`simulate.py` builds small Nano-BID-Amp datasets with known answers, for the WP2–WP4 tests and
the `-profile test` run. It needs Python ≥ 3.8 and `pyyaml` (`requirements.txt`).

```bash
python tests/simulate/simulate.py --design tests/simulate/design.yaml --outdir tests/data/synthetic
pytest tests/simulate            # MINIMAP2=/path/to/minimap2 also runs test_tools.py
```

The same design and seed give byte-identical output (`tests/data/synthetic/MANIFEST.md5`;
gzip headers carry no timestamp). `test_committed_data_is_current` fails if the committed data
drifts from the design, so edit `design.yaml` and regenerate rather than editing the data.

## Experiments

| Directory | Library | Samples | Purpose |
|---|---|---|---|
| `endogenous/` | endogenous, ONT adapters, UMI (defaults) | WT/KD × input/BS × 2 reps | default endogenous path |
| `mpra/` | MPRA, ONT + pool adapters, UMI (defaults) | WT/KD × input/BS × 2 reps | default MPRA path, near-duplicate oligos |
| `legacy_pool_noumi/` | MPRA, `orientation_adapters=pool`, `umi=false` (D29) | WT × input/BS × 2 reps | opt-in legacy mode (like paper run 20241114) |
| `legacy_3prime/` | endogenous, `ont_adapter_mode=three_prime_only` (D31) | WT × input/BS × 2 reps | opt-in legacy mode (like paper run 20250418); sense reads only |

Each directory has `samplesheet.csv` (§5.1), `analyses.yaml` (§5.4), `params.yaml` (run from the
repo root with `-params-file`), `fastq/` and `truth/`. References and BEDs are in `reference/`:
`<ref>.fa`, `<ref>.sites.bed` (single sites) and `<ref>.amplicons.bed` (whole targets, D11), all
standard 0-based (`bed_coordinates=bed0`).

## Design

- **Targets.** 5 endogenous amplicons (220–300 nt) and 20 MPRA oligos (130 nt) of random
  sequence without long homopolymers or adapter matches. `NEARDUP_B` is `NEARDUP_A` plus one
  base, so reads from A align equally well to both (MAPQ 0).
- **Sites.** Every site is a T flanked by non-T bases, and inserted bases differ from both
  neighbours, so each indel has a single alignment. Site types (`design.yaml: site_types`):
  `PUS7DEP`, `INDEP`, `UNMOD`, `INPUT_ZERO` (D25), `INS` (insertion after the T, D9), `LOWCOV`
  (20 vs 21 molecules: the `totalReads > 20` boundary) and `NOISY`.
- **Exact rates.** For each sample and site, `round(p × n)` molecules carry the deletion
  (insertion, substitution), so expected delrates have no sampling noise. Rep 2 uses a 0.85
  multiplier unless a per-rep list is given.
- **Molecules and duplicates.** Each molecule has a 10-nt UMI (pairwise Hamming ≥ 3 within a
  sample and target) and 1–3 PCR copies; half of the 3-copy molecules have one copy with a
  1-mismatch UMI, which UMICollapse must merge. All copies share the molecule's site outcomes, so
  dedup counts do not depend on which copy is kept.
- **Orientation.** Reads are sense or antisense at random (sense only in `legacy_3prime`).
- **Junk reads** per sample: `no_adapter`, `one_adapter_5p`, `one_adapter_3p` (dropped under
  D4), `too_short`, `short_pool` (passes the ONT trim, fails the pool trim length), `adapter_dimer`,
  `unmappable` (adapters around random sequence) and `chimera` (a sense and an antisense molecule
  in one read; counted in `n_reads_in_both_orientations`). Multimapper reads come from
  `NEARDUP_A`.
- **Error model.** Off (`error_model: none`): reads have no sequencing errors and quality 40.
  A Badread mode is reserved but not implemented.

## Truth tables (`<experiment>/truth/`)

| File | Content |
|---|---|
| `funnel.tsv` | `sample_id, step, records` for exactly the steps the mode runs (§5.2). `mapped` is given as the number of reads entering the aligner; minimap2 may add secondary records for the near-duplicate pair. |
| `both_orientations.tsv` | Expected `n_reads_in_both_orientations` (chimeras). |
| `site_counts.tsv` | Expected per-site counts in the §5.3 column order, after dedup (`level=dedup`) or per read when `umi=false` (`level=reads`), with `passes_min_coverage`. Rows with `passes_min_coverage=0` are the ones the pipeline drops. |
| `amplicon_counts.tsv` | The same for every T of every target (whole-target BED). |
| `site_counts_prededup.tsv` | Per-read counts before dedup (UMI experiments). |
| `categories.tsv` | Expected call per analysis and site; `clear_cut=1` rows are by design (`Modified`/`Unmodified` for `WT_mod`, `significant`/`not_significant` for `PUS7_dep_WT_v_KD`). |
| `reads.tsv.gz` | Every read: class, orientation, molecule(s), UMI(s), whether the UMI is the 1-mismatch copy, and fate (`mapped`, `mapq0`, `unmapped`, `dropped_trim`). |

Trimming truth comes from emulating the cutadapt passes on the actual reads (adapters are exact,
so exact matching is what cutadapt does with `-O 15`/`-O 10` and both adapters required).
Mapping truth is by construction and checked with minimap2 in `test_tools.py` (default presets,
D5): every valid record maps uniquely at position 1 with MAPQ ≥ 30 and exactly the designed
indels; near-duplicate records get MAPQ < 30; random inserts are unmapped.
