# Progress

Plan: [plan.md](plan.md) · Decisions: [docs/DECISIONS.md](docs/DECISIONS.md) · Change Register: [CHANGE_REGISTER.md](CHANGE_REGISTER.md)

| WP | Branch | Status | Notes |
|---|---|---|---|
| WP0 Scaffold and contracts | `wp/0-scaffold` | **Ready for G0 review** | nf-core template (tools 4.1.0), all §3 params incl. D29/D31, schemas, contracts, stubbed processes; `-profile test -stub` passes on Nextflow 25.04.7 in all three preprocessing modes and counting-only; validation errors checked; `nf-core pipelines lint` 0 failures (documented ignores in `.nf-core.yml`); CI in `.github/workflows/ci.yml` (nf-test not run locally: not installed on Sherlock) |
| WP1 Simulator | — | Not started | Can start after G0 |
| WP2 Preprocessing | — | Not started | Stub modules in `modules/local/`; D29/D31 branches wired in `subworkflows/local/preprocess.nf` |
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
2. **`bed_coordinates` default** (R-15, G1-e): currently required with no default.
3. **Golden endogenous across two runs.** Preprocessing parameters are per run, and the golden endogenous samples need different trimming per run (20250418 `three_prime_only`, 20251022 `linked`), so `test_golden` has separate `endogenous_20250418` and `endogenous_20251022` experiments. Site calling across both runs needs either a counts-level merge step (e.g. `--counts_merged` input for calling only) or per-sample preprocessing parameters in the sample sheet. Decide at G0.
4. **Containers.** Stub modules use a placeholder `ubuntu:22.04` image; WP2–WP4 pin real containers.
5. **Boolean CLI flags on Nextflow edge.** Nextflow 26.09.2-edge rejects `--umi false` given on the command line ("Value is [string] but should be [boolean]"); 25.04.7 accepts it. CI runs `latest-everything` as non-blocking. If stable releases keep this behaviour, document `-params-file` (or `--umi=false` alternatives) for boolean options.
6. Endogenous random effects in `analyses_example.yaml` (`(1|vector)` vs `(1|celltype)`; plan §10).
