# Progress

Plan: [plan.md](plan.md) · Decisions: [docs/DECISIONS.md](docs/DECISIONS.md) · Change Register: [CHANGE_REGISTER.md](CHANGE_REGISTER.md)

| WP | Branch | Status | Notes |
|---|---|---|---|
| WP0 Scaffold and contracts | `wp/0-scaffold` | **Ready for G0 review** | nf-core template (tools 4.1.0), all §3 params incl. D29/D31, schemas, contracts, stubbed processes; `-profile test -stub` passes on Nextflow 25.04.7 in all three preprocessing modes and counting-only; validation errors checked; `nf-core pipelines lint` 0 failures (documented ignores in `.nf-core.yml`); CI in `.github/workflows/ci.yml` (nf-test not run locally: not installed on Sherlock) |
| WP1 Simulator | — | Not started | Can start after G0 |
| WP2 Preprocessing | — | Not started | Stub modules in `modules/local/`; D29/D31 branches wired in `subworkflows/local/preprocess.nf` |
| WP3 Counting | `wp/3-count` (local, not pushed) | **Implemented; local tests pass** | R package `R/nanomodamp` (`count_sites`, `merge_counts`, `merge_count_tables` (D32), `background_qc`), CLIs `bin/nma_count.R`, `nma_merge.R`, `nma_bgqc.R`, modules wired (stubs kept). 23 testthat tests pass on R 4.3.2 / Bioc 3.18, incl. exact equivalence with the legacy `bam_counts_fast.R` on a fixture BAM. On golden V2 BAMs with the 0-based `targets.bed` (`bed0`), counts are identical to `expected/legacy_rerun` for in vitro 4/4, in cellulo 32/32, endogenous 20251022 16/16. See "WP3 — remaining". |
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
4. **Containers.** Stub modules use a placeholder `ubuntu:22.04` image; WP2–WP4 pin real containers.
5. **Boolean CLI flags on Nextflow edge.** Nextflow 26.09.2-edge rejects `--umi false` given on the command line ("Value is [string] but should be [boolean]"); 25.04.7 accepts it. CI runs `latest-everything` as non-blocking. If stable releases keep this behaviour, document `-params-file` (or `--umi=false` alternatives) for boolean options.
6. ~~Endogenous random effects~~ — decided 2026-10-09 (D33): `(1|rep)` + `(1|celltype)` when both cell types are combined.

## WP3 — remaining (2026-10-09)
- **Container:** `containers/nanomodamp-r/Dockerfile` written but not built or pushed; modules still use the `ubuntu:22.04` placeholder, so `-profile docker` real runs of COUNT/MERGE/BGQC will fail until the image is built and pinned (needs a registry decision from Becca).
- **CI:** the R testthat job (`devtools::test`) has not run yet (branch not pushed). Locally the suite was run with `testthat::test_dir()` on an installed package.
- **nf-test:** no module-level nf-test for COUNT/MERGE yet (nf-test not installed on Sherlock); a real `--input_counts` run through Nextflow was checked by hand (merge correct; MULTIQC failed only because multiqc isn't installed locally).
- **Golden 20250418:** not compared with `count_sites()` here — its legacy counts came from the older `bam_counts.R` on pre-dedup BAMs; WP7 covers it.
- **Note for WP2/WP7:** a stale `.fai` next to `--fasta` is used as is (the pipeline only builds a missing one); the paper's pool1 `.fai` is stale. Consider always rebuilding or validating the `.fai`.
