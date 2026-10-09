# Progress

Plan: [plan.md](plan.md) · Decisions: [docs/DECISIONS.md](docs/DECISIONS.md) · Change Register: [CHANGE_REGISTER.md](CHANGE_REGISTER.md)

| WP | Branch | Status | Notes |
|---|---|---|---|
| WP0 Scaffold and contracts | `wp/0-scaffold` | **Ready for G0 review** | nf-core template (tools 4.1.0), all §3 params incl. D29/D31, schemas, contracts, stubbed processes; `-profile test -stub` passes on Nextflow 25.04.7 in all three preprocessing modes and counting-only; validation errors checked; `nf-core pipelines lint` 0 failures (documented ignores in `.nf-core.yml`); CI in `.github/workflows/ci.yml` (nf-test not run locally: not installed on Sherlock) |
| WP1 Simulator | — | Not started | Can start after G0 |
| WP2 Preprocessing | — | Not started | Stub modules in `modules/local/`; D29/D31 branches wired in `subworkflows/local/preprocess.nf` |
| WP3 Counting | — | Not started | |
| WP4 Site calling | `wp/4-calling` (local, not pushed) | **Implemented; testthat 1–11 pass locally** | `R/nanomodamp` calling/sitesets/plots/run; `bin/nma_call.R`, `bin/nma_sitesets.R`; CALL_SITES / SITE_SETS modules call them (stub run passes); `assets/analyses_example.yaml` + `assets/analyses/{invitro,incellulo,endogenous}.yaml`. Test 11: exact equivalence with the legacy scripts rerun on the shipped Figure 3 tables (see below). New discrepancy **R-30 (PROPOSED)**: legacy BH is a no-op |
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

## WP4 — status and remaining (2026-10-09)

**Test 11 (legacy equivalence)** — references are `modification_analysis.R` / `incell_analysis.R`
rerun unchanged (R 4.3.2) on the shipped Figure 3 count tables
(`R/nanomodamp/tests/testthat/fixtures/legacy_equivalence/`):
- In vitro, all 770 sites: identical order, categories and TOST labels; p-values max relative
  difference 5.9e-15; deltas/eqbounds ≤ 4e-16. Shipped paper table: same categories (p within 3e-4
  relative: it came from slightly different script/package versions).
- In cellulo WT_mod (HepG2, 293T, Both) and all six PUS7_dep runs + union: identical on the
  default deterministic subset (60 sites) — see the full-run line below. One p-value
  (PFKP_chr10_3112271, 6.5e-315) is a readr serialization artifact of subnormal doubles.
- `bin/nma_call.R` on the in vitro table reproduces the legacy table exactly.

**Open for Becca / Orchestrator**
1. **R-30 (PROPOSED):** both legacy scripts run `p.adjust(..., "BH")` inside a `mutate()` on a
   tibble still grouped by `(chr, pos)`, so no multiple-testing correction is applied (shipped
   tables: `p.value == p.adjust_diff == p.value.BH`). Ported as `p_adjust = "legacy"` (default,
   reproduces the paper); `p_adjust = "BH"` implements plan §6.3/§6.4 as written. Not yet exposed
   in the analyses YAML / `schema_analyses.json` (contract change).
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
