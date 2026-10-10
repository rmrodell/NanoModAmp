# Progress

Plan: [plan.md](plan.md) · Decisions: [docs/DECISIONS.md](docs/DECISIONS.md) · Change Register: [CHANGE_REGISTER.md](CHANGE_REGISTER.md)

| WP | Branch | Status | Notes |
|---|---|---|---|
| WP0 Scaffold and contracts | `wp/0-scaffold` | **Done — G0 approved 2026-10-09** | nf-core template, schemas, contracts, CI; merged into `dev` |
| WP1 Simulator | `wp/1-simulator` | **Merged into `dev`** | `tests/simulate/`, `tests/data/synthetic/` (1.2 MB); pytest 41 passed on Sherlock |
| WP2 Preprocessing | `wp/2-preprocess` | **Merged into `dev`; CI pending** | Real modules, all modes verified on Sherlock (Apptainer); nf-test specs first run in CI. See "WP2 — remaining" |
| WP3 Counting | `wp/3-count` | **Merged into `dev`** | `R/nanomodamp` counting/merge/bg-QC, R-31 fresh `.fai`; 24 testthat tests; golden counts identical to `legacy_rerun`. See "WP3 — remaining" |
| WP4 Site calling | `wp/4-calling` | **Merged into `dev`** | Treatment/factor tests, site sets, plots; test 11 exact vs paper scripts (in vitro 770, in cellulo 760 sites); R-30 BH default. See "WP4 — remaining" |
| WP5 Golden package | `wp/5-golden-package` | **Merged into `dev`** | `tests/data/golden/` (22.3 MB); L1, L2, V1, V2, Repro pass; build log `dev/golden/BUILD_LOG.md` |
| WP6 Documentation | — | Started in WP0 | `docs/usage.md` legacy-library section, `docs/methods.md` stub |
| WP7 Integration | — | Blocked on WP1–WP6 | Compares against `legacy_rerun`, `published_rerun`, `paper_reference` (D30) |
| WP8 Release | — | Blocked | |

## Blockers

- G0 approved 2026-10-09 (Becca); contracts are frozen — changes go through `docs/contracts/CHANGELOG.md`.

## Open questions (raised for G0; G0 approved 2026-10-09)

1. **Nextflow minimum version.** The nf-core 4.1.0 template requires Nextflow ≥ 25.10.4; Sherlock's newest module is 25.04.7. WP0 relaxed the pin to `>=25.04.7` and verified the stub runs there. Keep 25.04.7 as the minimum, or require 25.10 (then Sherlock users need their own Nextflow install)?
2. ~~`bed_coordinates` default~~ — **decided 2026-10-06 (G1-e, D27/R-15):** default `bed0`; `one_based_start` for legacy BEDs; start = end rows under `bed0` fail validation.
3. ~~Golden endogenous across two runs~~ — **decided 2026-10-06 (D32):** merge at the counts level with `--input_counts`; implemented as a stub (`MERGE_COUNT_TABLES`, WP3 implements it; WP2/WP3 test 15).
4. **Containers.** Stub modules use a placeholder `ubuntu:22.04` image; WP2–WP4 pin real containers. WP2 done: preprocessing modules pin biocontainers / Seqera community images (see each module's `container`).
5. **Boolean CLI flags on Nextflow edge.** Nextflow 26.09.2-edge rejects `--umi false` given on the command line ("Value is [string] but should be [boolean]"); 25.04.7 accepts it. CI runs `latest-everything` as non-blocking. If stable releases keep this behaviour, document `-params-file` (or `--umi=false` alternatives) for boolean options.
6. ~~Endogenous random effects~~ — decided 2026-10-09 (D33): `(1|rep)` + `(1|celltype)` when both cell types are combined.

## WP3 — remaining (2026-10-09)
- **Container:** `containers/nanomodamp-r/Dockerfile` written but not built or pushed; modules still use the `ubuntu:22.04` placeholder, so `-profile docker` real runs of COUNT/MERGE/BGQC will fail until the image is built and pinned (needs a registry decision from Becca).
- **CI:** the R testthat job (`devtools::test`) has not run yet (branch not pushed). Locally the suite was run with `testthat::test_dir()` on an installed package.
- **nf-test:** no module-level nf-test for COUNT/MERGE yet (nf-test not installed on Sherlock); a real `--input_counts` run through Nextflow was checked by hand (merge correct; MULTIQC failed only because multiqc isn't installed locally).
- **Golden 20250418:** not compared with `count_sites()` here — its legacy counts came from the older `bam_counts.R` on pre-dedup BAMs; WP7 covers it.
- ~~Stale `.fai` used as is~~ — fixed (R-31, 2026-10-09): SAMTOOLS_FAIDX always rebuilds the index on a copy of `--fasta`; `check_fai()` in counting rejects a mismatched index.
## WP2 — remaining (2026-10-09)

- **Run the nf-test specs in CI** (`subworkflows/local/preprocess/tests/*.nf.test`, `modules/local/cat_fastq/tests/main.nf.test`); they were checked locally only through the equivalent harness (`tests/fixtures/preprocess/harness.nf` + `check_outputs.py`, all checks pass in all four modes) and direct tool runs. Requires pushing `wp/2-preprocess` (Becca's approval).
- **Test 11** (funnel equals the WP1 truth): currently checked against the fixture truth in `tests/fixtures/preprocess/expected.json`; repeat on the WP1 synthetic data once WP1 merges.
- **Test 12** (`-profile test` end to end, both library types, real tools): needs WP3 counting and WP1 data; stub runs of `-profile test -stub` pass in all three modes.
- **Test 9 (missing file)** is enforced by the sample-sheet schema at pipeline level (WP0 test); CAT_FASTQ also fails clearly on empty files and on directories without `*.fastq.gz`.
- **Observation for WP7 (R-04):** on the fixtures, reads carrying only one ONT adapter are discarded both by cutadapt 5.2 with `;required` and by legacy cutadapt 1.18 linked `-g` with `--discard-untrimmed`; the R-04 difference will have to be measured on real data (golden package).
- **Lint:** `nextflow lint` passes on all new files; `nf-core pipelines lint` did not finish locally (stalled after generating container configs, 15 min timeout) — CI runs it.
- **Fixture observation:** the chimera fixture (sense + antisense construct in one read) puts 1 of 9 merged reads in both orientations, so MERGE_ORIENT logs the > 0.1% warning there by design.
## WP4 — status and remaining (2026-10-09)
**Test 11 (legacy equivalence)** — references are `modification_analysis.R` / `incell_analysis.R`
rerun unchanged (R 4.3.2) on the shipped Figure 3 count tables
(`R/nanomodamp/tests/testthat/fixtures/legacy_equivalence/`):
- In vitro, all 770 sites: identical order, categories and TOST labels; p-values max relative
  difference 5.9e-15; deltas/eqbounds ≤ 4e-16. Shipped paper table: same categories (p within 3e-4
  relative: it came from slightly different script/package versions).
- In cellulo WT_mod (HepG2, 293T, Both) and all six PUS7_dep runs + union: identical to the
  legacy rerun on **all 760 sites** (full run, `NANOMODAMP_FULL_EQUIVALENCE=true`, 34 min on one
  shared CPU) and on the default 60-site subset (≈3 min). One p-value (PFKP_chr10_3112271,
  6.5e-315) is a readr serialization artifact of subnormal doubles.
- Paper (shipped) vs legacy rerun: PUS7 union identical (184 sites); WT_mod_Both categories
  identical except RPL22_chr1_6186768 (paper Inconclusive, rerun Unmodified; p 1.4e-5 vs 3.0e-22),
  i.e. drift between the original run and today's packages, not a porting difference. Documented
  in the fixture README and allowed explicitly in the test.
- `bin/nma_call.R` on the in vitro table reproduces the legacy table exactly.
**Open for Becca / Orchestrator**
1. ~~R-30~~ — decided by Becca 2026-10-09: default `p_adjust: BH` (schema key added, contract CHANGELOG updated); `legacy` reproduces the paper (used by test 11).
2. Treatment table appends `p.value`, `is_equivalent`, `all_below_thresh` (contract §5.5 lists
   `p.value`, `is_equivalent`, `model_status` as allowed appends); model failures stay in
   `equivalence_status` as in the legacy code. Orchestrator: accept `all_below_thresh` or drop it.
3. In-cellulo legacy script labels equivalence `Equivalent`/`Not equivalent`; the port uses the
   detailed `modification_analysis.R` labels everywhere (plan §6.3); categories are unaffected.
**WP4 — remaining**
- Container: CALL_SITES / SITE_SETS still use the placeholder image; needs the `nanomodamp-r`
  image (R ≥ 4.3 + nanomodamp + blme/lme4/car/furrr/ggplot2/forcats/yaml/argparse), shared with WP3.
- CI: run testthat with `NANOMODAMP_FULL_EQUIVALENCE=true` in a dedicated job (~1 h on one CPU;
  parallel workers do not change results).
- Merge with WP3: both branches add `R/nanomodamp/DESCRIPTION` and `NAMESPACE`; resolve by taking
  the union of Imports and export()/import lines.
- nf-test for CALL_SITES / SITE_SETS with real (non-stub) runs once the container exists.
