# NanoModAmp golden test package

A small public subset of the paper data (*RNA sequence, structure, and cell type specific
features drive pseudouridylation by PUS7*, Figures 2 and 3). Run it to check that an
installation reproduces known results. CI uses it as a regression test.

Built by `dev/golden/assemble_package.py` (plan: `docs/plans/golden_test_package_plan.md`;
decisions: `docs/DECISIONS.md`). `provenance.json` records inputs (basenames + md5), tool
versions, seed and per-sample read counts; `MANIFEST.md5` lists every file.

## Contents

| Directory | Data | Samples | Targets |
|---|---|---|---|
| `endogenous/` | Figure 2, runs 20250418 and 20251022 | 32 (16 per run) | RHBDD2, HDAC6 (20250418); STIM1 (20251022) |
| `mpra_invitro/` | Figure 3 in vitro, May run | 4 | 6 pool1 oligos (`pool1_chr.txt`) |
| `mpra_incell/` | Figure 3 in cellulo, run 20241114 | 32 | same 6 oligos |

Each directory has `fastq/` (raw reads, adapters intact, original order), `reference/` (the
**full** reference, bgzip + index), `targets.bed` (paper single-site BED filtered to the targets,
**converted to standard 0-based**: start = site − 1, end = site; run with `bed_coordinates=bed0`), sample sheet(s), `params*.yaml`, `analyses.yaml` and `expected/`.

## How the reads were chosen (G1-b, G1-d)

All reads with any alignment on a target in the legacy rerun's pre-dedup BAM, with their PCR
duplicates, then at most **250 reads per sample and target** (seeded, nested) plus **150 random
background reads per sample** (seeded) so the adapter, length and mapping filters are exercised.
Targets below the cap are kept whole. Seed 20261004.

## Running it

    nextflow run . -profile test_golden,<docker|singularity> --golden_experiment <name> --outdir <dir>

`<name>` is `endogenous_20250418`, `endogenous_20251022`, `mpra_invitro` or `mpra_incell`.
Endogenous site calling runs on the counts-level merge of both runs (D32):
`--input_counts <dir_20250418>/counts/counts_merged.tsv,<dir_20251022>/counts/counts_merged.tsv
--analyses tests/data/golden/endogenous/analyses.yaml` (see `endogenous/params_merged.yaml`).

### Non-default modes — these samples only

Pipeline defaults (ONT trim with both adapters, UMI extraction, deduplication) apply to all
data, including in vitro and in cellulo. Two paper libraries need opt-in settings, set only for
them in `conf/test_golden.config` and the `params*.yaml` here:

- **Run 20250418 (endogenous):** reads lack the 5′ ONT adapter → `ont_adapter_mode=three_prime_only` (D31, R-29).
- **Run 20241114 (in cellulo):** library built without ONT adapters or UMIs → `orientation_adapters=pool`, `umi=false` (D29, R-26).

### Sample-sheet conventions

- Endogenous vectors are normalised: P36/pLKO → WT, P97/shPUS7 → KD; reps are run-unique
  (`0418r1`, `1022r2`, …) because batch pairing holds only within a run. `endogenous/name_map.tsv`
  maps package names to the original sample names.
- In vitro noPUS samples are `treat=input` (D14), noted in `notes`.

## Expected outputs: what is binding

| Directory | Source | Status |
|---|---|---|
| `expected/legacy_rerun/` | Scripts that actually produced the paper numbers (`dev/legacy/`), rerun on this package (V2) | **Binding for counts.** Endogenous calls are reference only (legacy used `glmer`; the pipeline uses `bglmer`, D13). |
| `expected/published_rerun/` | Published Figure 2/3 scripts (`dev/published/`), rerun on this package (D30) | Reference; isolates implementation differences in WP7 |
| `expected/legacy_full_selected/` | Legacy counts on the selection before the read cap (V1) | Reference; V1 matched the paper counts |
| `expected/paper_reference/` | Paper count tables and shipped site-call tables, restricted to these samples and targets | Reference only: paper calls used all samples and sites |

Notes:
- The paper processing differs from the published scripts (R-25): minimap2 `-a -k5`, no
  primary/MAPQ filter, `umi_tools dedup --method directional`. The legacy harness adds a fixed
  `--random-seed` (G1-c), so `legacy_rerun` is deterministic but its dedup representatives
  differ from the paper's.
- Paper counts for 20250418 were made on **pre-dedup** BAMs; 20251022 and in vitro on dedup BAMs;
  in cellulo has no dedup (R-27). V1 was exact for non-deduplicated counts and within tolerance
  for deduplicated ones (G1-c).
- `endoG1G2.fasta` is used for all endogenous runs (F2); target sequences are identical to the
  per-run references the paper used.
- BED convention (R-15, gate G1-e): the pipeline default is `bed_coordinates=bed0`. The paper's
  delpos BEDs have start = end = the 1-based site (the legacy counter adds no +1); the package
  `targets.bed` is converted to standard 0-based. The legacy and published reruns used the
  original start = end BED with the legacy code. Positions in all expected outputs are 1-based
  sites and are identical under both conventions.
- Run 20250418 paper counts also contain window-start rows from its whole-amplicon BED; these
  are a documented paper artifact (`dev/golden/R15_bed_convention.md`) and are not in the package.
- Known gaps (G1-b): no Inconclusive, input-delrate = 0 or Modified-not-PUS7-dependent MPRA
  site; cell-type-specific analyses are out of scope.
