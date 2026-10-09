# Progress

Plan: [plan.md](plan.md) · Decisions: [docs/DECISIONS.md](docs/DECISIONS.md) · Change Register: [CHANGE_REGISTER.md](CHANGE_REGISTER.md)

| WP | Branch | Status | Notes |
|---|---|---|---|
| WP0 Scaffold and contracts | `wp/0-scaffold` | **Ready for G0 review** | nf-core template (tools 4.1.0), all §3 params incl. D29/D31, schemas, contracts, stubbed processes; `-profile test -stub` passes on Nextflow 25.04.7 in all three preprocessing modes and counting-only; validation errors checked; `nf-core pipelines lint` 0 failures (documented ignores in `.nf-core.yml`); CI in `.github/workflows/ci.yml` (nf-test not run locally: not installed on Sherlock) |
| WP1 Simulator | — | Not started | Can start after G0 |
| WP2 Preprocessing | `wp/2-preprocess` (local, not pushed) | **Implemented; CI not yet run** | Real modules (nf-core cutadapt 5.2, seqtk 1.4, umi_tools, minimap2 + samtools 1.24, umicollapse 1.1.0; local CAT_FASTQ, MERGE_ORIENT, READ_FUNNEL); all three modes (D29, D31) verified on Sherlock with Apptainer on hand-made fixtures; nf-test specs for tests 1–11, 13, 14 written (not run: nf-test not installed on Sherlock). See "WP2 — remaining" |
| WP3 Counting | — | Not started | |
| WP4 Site calling | — | Not started | Test 11 (legacy equivalence on shipped Figure 3 tables) is independent of WP5 |
| WP5 Golden package | `wp/5-golden-package` (local, not pushed) | In progress | Harnesses, selection, L1/V1/V2 and published rerun running on Sherlock; see `docs/plans/golden_test_package_plan.md` |
| WP6 Documentation | — | Started in WP0 | `docs/usage.md` legacy-library section, `docs/methods.md` stub |
| WP7 Integration | — | Blocked on WP1–WP6 | Compares against `legacy_rerun`, `published_rerun`, `paper_reference` (D30) |
| WP8 Release | — | Blocked | |

## Blockers

- None for WP0. WP1–WP6 wait for G0 (contracts frozen).

## Open questions (for G0)

1. **Nextflow minimum version.** The nf-core 4.1.0 template requires Nextflow ≥ 25.10.4; Sherlock's newest module is 25.04.7. WP0 relaxed the pin to `>=25.04.7` and verified the stub runs there. Keep 25.04.7 as the minimum, or require 25.10 (then Sherlock users need their own Nextflow install)?
2. ~~`bed_coordinates` default~~ — **decided 2026-10-06 (G1-e, D27/R-15):** default `bed0`; `one_based_start` for legacy BEDs; start = end rows under `bed0` fail validation.
3. ~~Golden endogenous across two runs~~ — **decided 2026-10-06 (D32):** merge at the counts level with `--input_counts`; implemented as a stub (`MERGE_COUNT_TABLES`, WP3 implements it; WP2/WP3 test 15).
4. **Containers.** Stub modules use a placeholder `ubuntu:22.04` image; WP2–WP4 pin real containers. WP2 done: preprocessing modules pin biocontainers / Seqera community images (see each module's `container`).
5. **Boolean CLI flags on Nextflow edge.** Nextflow 26.09.2-edge rejects `--umi false` given on the command line ("Value is [string] but should be [boolean]"); 25.04.7 accepts it. CI runs `latest-everything` as non-blocking. If stable releases keep this behaviour, document `-params-file` (or `--umi=false` alternatives) for boolean options.
6. ~~Endogenous random effects~~ — decided 2026-10-09 (D33): `(1|rep)` + `(1|celltype)` when both cell types are combined.

## WP2 — remaining (2026-10-09)

- **Run the nf-test specs in CI** (`subworkflows/local/preprocess/tests/*.nf.test`, `modules/local/cat_fastq/tests/main.nf.test`); they were checked locally only through the equivalent harness (`tests/fixtures/preprocess/harness.nf` + `check_outputs.py`, all checks pass in all four modes) and direct tool runs. Requires pushing `wp/2-preprocess` (Becca's approval).
- **Test 11** (funnel equals the WP1 truth): currently checked against the fixture truth in `tests/fixtures/preprocess/expected.json`; repeat on the WP1 synthetic data once WP1 merges.
- **Test 12** (`-profile test` end to end, both library types, real tools): needs WP3 counting and WP1 data; stub runs of `-profile test -stub` pass in all three modes.
- **Test 9 (missing file)** is enforced by the sample-sheet schema at pipeline level (WP0 test); CAT_FASTQ also fails clearly on empty files and on directories without `*.fastq.gz`.
- **Observation for WP7 (R-04):** on the fixtures, reads carrying only one ONT adapter are discarded both by cutadapt 5.2 with `;required` and by legacy cutadapt 1.18 linked `-g` with `--discard-untrimmed`; the R-04 difference will have to be measured on real data (golden package).
- **Fixture observation:** the chimera fixture (sense + antisense construct in one read) puts 1 of 9 merged reads in both orientations, so MERGE_ORIENT logs the > 0.1% warning there by design.
