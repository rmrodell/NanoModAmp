# Change Register

Seeded from `plan.md` §9 (WP0, 2026-10-06). Every difference from the paper's results must be explained by an entry here (D20). New discrepancies are added as `PROPOSED` and raised at the next gate (plan §0 rule 1).

Status values: `APPROVED` (implement), `KEEP` (preserve legacy behavior and document it), `SUPERSEDED` (not applicable in the new design), `OPEN` (decide at a gate), `PROPOSED` (raised, awaiting a gate).

| ID | Legacy behavior (file) | Resolution | Status | Expected effect on results |
|---|---|---|---|---|
| R-01 | `sample_name.R` missing (`BIDdetect.sh`) | Sample sheet metadata replaces it; reconstructed for the legacy harness only (D26) | SUPERSEDED | None |
| R-02 | BH applied over the `(Intercept)` and `treat` rows (`analysis_endo.R`) | Replaced by the Figure 3 treatment analysis (D13) | SUPERSEDED | Endogenous p.adj/calls may change; report in WP7 |
| R-03 | `--cores` does nothing: `purrr::map` used instead of `future_map` (`modification_analysis.R`) | Use seeded `future_map` | APPROVED | None (performance only) |
| R-04 | Linked adapters with `-g` are non-anchored and may accept reads with one flank; cutadapt 1.18 | Require both flanks (D4); modern pinned cutadapt | APPROVED | Fewer reads; small count differences. Quantify in WP7 |
| R-05 | Master counts file appended on rerun (`BIDdetect.sh`) | Per-sample outputs plus a name-based merge | SUPERSEDED | None |
| R-06 | Columns combined by position under a hard-coded header (`BIDdetect.sh`) | Join by name (§5.3) | APPROVED | None |
| R-07 | Missing `/` in input path (`analysis_endo.R`) | N/A | SUPERSEDED | None |
| R-08 | Raw `.fastq.gz` count reported as `N/A` (`track_metrics`) | Count gz correctly | APPROVED | Metrics only |
| R-09 | `wc -l` sample map; blank lines, missing trailing newline | Schema-validated sample sheet | SUPERSEDED | None |
| R-10 | `REF_FA` unbound under `set -u`; help/default time mismatch (`submit_file_prep.sh`) | Replaced by Nextflow | SUPERSEDED | None |
| R-11 | `calculate_delta_delrate` returns a string when there are no sites | N/A | SUPERSEDED | None |
| R-12 | Per-region errors swallowed with warnings (`bam_counts_fast.R`) | Fail the task unless `allow_region_failures`; always write `failed_regions.tsv` | APPROVED | None if no failures |
| R-13 | Coverage `> 20` at counting vs `>= 20` in the sweep wording | Keep `> 20` (D10); document | KEEP | None |
| R-14 | Cosmetic: step numbering, "umi_tools dedup" label, "cleanup disabled" message, README says MAPQ > 30 | N/A | SUPERSEDED | None |
| R-15 | BED start used without the 0→1 conversion (commented out); the paper's single-site BEDs were written as start = end = 1-based site to match | `bed_coordinates` (default `bed0`; `one_based_start` for legacy BEDs). Golden `targets.bed` converted to 0-based. The 20250418 window-BED rows at the window start (e.g. RHBDD2:285) were not intended sites: a known paper artifact, not in the golden package | APPROVED (G1-e, 2026-10-06) | None for the golden targets; legacy BEDs read with `bed0` would shift by one |
| R-16 | `max_depth=200000` may truncate silently | Keep the default; warn when reached | APPROVED | None unless saturated |
| R-17 | `all_below_thresh` `na.rm` differs between scripts | Use `na.rm=TRUE` (the in-cellulo version) | APPROVED | None expected (delrate is never NA) |
| R-18 | Insertions included in `totalReads` | Keep (D9); document | KEEP | None |
| R-19 | Equivalence bound infinite when input mean = 0 | Keep (D25); document | KEEP | None |
| R-20 | Some plots only as PDF (all-sites, labeled heatmap) | Always PDF + PNG (D3) | APPROVED | Outputs only |
| R-21 | Endogenous analysis used `lme4::glmer` with vector random effect | `bglmer` Figure 3 method for all (D13); endogenous random effects set via config | APPROVED | Endogenous calls may change; report in WP7 |
| R-22 | Endogenous Δ filter before testing and BH | Figure 3 order (D16) | APPROVED | Endogenous calls may change; report in WP7 |
| R-23 | In-vitro `vector=="noPUS"` recoded to input | Explicit `treat` in the sample sheet (D14) | APPROVED | None if the golden sheet labels them correctly |
| R-24 | Hard-coded tool paths, unversioned modules, runtime installs, `detectCores()` | Containers, pinned versions, `task.cpus` | APPROVED | None |
| R-25 | Paper data were processed with earlier scripts: minimap2 `-a -k5`, no primary/MAPQ filter, `umi_tools dedup --method directional` (not the published `-ax sr`/`splice -uf`, `-F 2304 -q 30`, UMICollapse) | Pipeline follows the published methods (D5, D7). Effect measured as `legacy_rerun` vs `published_rerun` on the golden package (D30) | APPROVED | Counts and calls differ from paper numbers; quantify in WP7 |
| R-26 | The paper's in-cellulo samples (run 20241114) came from a library without ONT adapters or UMIs; the paper used a separate script (pool-adapter orientation, no dedup) | Opt-in parameters `orientation_adapters=pool`, `umi=false` (D29), set only for those golden samples; defaults (ONT + UMI + dedup) apply to all other data, including new in-cellulo data | APPROVED | None vs the in-cellulo legacy design |
| R-27 | Endogenous 20250418 paper counts were made on **pre-dedup** BAMs; 20251022 on dedup BAMs | Pipeline always counts the final BAM; difference reported for 20250418 | APPROVED | 20250418 totalReads higher in the paper than after dedup |
| R-29 | Paper run 20250418 reads lack the 5′ ONT adapter; legacy trimmed the 3′ adapter only (`-a`, untrimmed kept) | Opt-in `ont_adapter_mode=three_prime_only` (D31) for those golden samples only; published rerun uses the same trim for them (documented variant) | APPROVED | Without it, 0 reads pass for 20250418 |
| R-30 | BH correction is a no-op: `p.adjust(p.value, "BH")` runs inside `mutate()` on a tibble still grouped by `(chr, pos)` (`modification_analysis.R::diff_equiv_analysis`, `incell_analysis.R::{diff_equiv_analysis, perform_anova_test}`), so `p.adjust_diff` / `p.value.BH` equal the raw p-values | Default BH; legacy option reproduces paper (`p_adjust: legacy`, WP4 test 11) | APPROVED (Becca, 2026-10-09) | Fewer Modified / PUS7-dependent calls than the paper; quantify in WP7 |
| R-28 | `umi_tools dedup` picks random representatives (no seed) | Superseded by UMICollapse (D7); legacy harness seeded (edit e) | SUPERSEDED | None |

---
