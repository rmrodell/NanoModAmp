# Legacy-equivalence fixtures (WP4 testthat case 11)

Source: `martinezlab/PUS7regulation2026` @ 4fd285b, Figure 3.

| File | Origin | md5 |
|---|---|---|
| `legacy_scripts/modification_analysis.R` | `Figure3/mpra_sites/modification_analysis.R`, byte-identical | 8abf253f160547b5a046a9ef0e04aba3 |
| `legacy_scripts/incell_analysis.R` | `Figure3/mpra_sites/incell_analysis.R`, byte-identical (CRLF line endings; `-text` in `.gitattributes`) | dd76753112f15a4746f1869b054df449 |
| `inputs/BIDdetect_data_invitro_delpos.txt.gz` | `Figure3/parameter_sweep/BIDdetect_data_invitro_delpos.txt` (gzip -9 -n) | source 6abb0ec74cbba068248e47c875f82316 |
| `inputs/BIDdetect_data_incell_delpos.txt.gz` | `Figure3/plots/BIDdetect_data_incell_delpos.txt` (gzip -9 -n) | source 9eacfaeb3b2dc4f33cc9eddf71f1e0db |
| `shipped/*.tsv` | `Figure3/plots/` paper outputs, unchanged | — |
| `legacy_rerun/invitro/` | `modification_analysis.R` rerun on Sherlock on the in-vitro input: `--prefix invitro_delpos --cores 1` | — |
| `legacy_rerun/incell/` | `incell_analysis.R` rerun on Sherlock on the in-cellulo input; only the user-configuration block changed (input/output paths, `num_cores <- 1`, `plot_all_sites <- FALSE`, `run_cell_type_analysis <- FALSE`) | — |

## How test 11 compares (since 2026-11-05)

GLMM fits (`blme::bglmer`) give slightly different estimates under different R / lme4 / blme
versions: on GitHub's newest R the p-values moved by 1e-4 to 1e-2 relative from the Sherlock
references and one in-vitro site changed category. So test 11 compares **like with like**: at test
time it runs the vendored legacy scripts (`Rscript` of the running R, same library) and the port on
the same input, and asserts p-values to 1e-8 relative and categories, labels, deltas and site sets
exactly. `incell_analysis.R` is configured by editing a temporary copy of its user-configuration
block (the same four edits as for the stored rerun); the test checks that each edit matched exactly
one line.

- **In vitro (11a):** all 770 sites (~1 min for the legacy script).
- **In cellulo (11b, 11c):** by default a deterministic subset of oligos (`incell_subset_chrs()` in
  the test: paper-vs-rerun drift sites, sites nearest the decision boundaries, the most and least
  significant called sites of every PUS7 run, plus 12 seeded random oligos). Per-site results don't
  depend on the subset: every model is fitted per site and the legacy "BH" is a per-site no-op
  (R-30). `NANOMODAMP_FULL_EQUIVALENCE=true` runs all 760 oligos (legacy ~50 min + port ~35 min,
  one core).

The stored Sherlock rerun (`legacy_rerun/`, R 4.3.2, lme4 1.1.35.5, blme 1.0.6, car 3.1.3, dplyr
1.1.4, 2026-10-09, `/scratch/users/rodell/NanoModAmp_build/wp4_legacy_ref/run_refs.sh`) and the
shipped paper tables stay as documented references: the test asserts against them only when the
running R has exactly those versions, otherwise it reports differences with `message()`.

## Findings while building these fixtures

- **R-30 (APPROVED 2026-10-09, default BH):** both legacy scripts apply `p.adjust(..., "BH")` inside
  a `mutate()` on a tibble still grouped by `(chr, pos)`, so the "adjusted" p-values equal the raw
  ones (visible in the shipped tables: `p.value == p.adjust_diff`, `p.value == p.value.BH`). The port
  reproduces this with `p_adjust = "legacy"`; test 11 uses it.
- readr writes subnormal doubles with a wrong exponent (e.g. 6.5187e-315 is written as
  `6.5187e-298`, PFKP_chr10_3112271 in `WT_mod_HepG2`); the test treats p < 1e-290 as equal.
- **Paper vs rerun drift (in cellulo WT_mod_Both):** rerunning the unchanged legacy script on
  Sherlock reproduces the shipped categories at 759 of 760 sites; RPL22_chr1_6186768 is Inconclusive
  in the paper and Unmodified in the rerun (p 1.4e-5 vs 3.0e-22). 23 of 341 tested p-values differ
  by more than 1e-3 relative. The PUS7 union is identical (184 sites). The port equals the rerun
  exactly, so the drift comes from the original run's package versions / optimizer, not the port.
