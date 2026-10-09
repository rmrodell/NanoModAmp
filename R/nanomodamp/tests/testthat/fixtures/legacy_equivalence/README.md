# Legacy-equivalence fixtures (WP4 testthat case 11)

Source: `martinezlab/PUS7regulation2026` @ 4fd285b, Figure 3.

| File | Origin | md5 of the source |
|---|---|---|
| `inputs/BIDdetect_data_invitro_delpos.txt.gz` | `Figure3/parameter_sweep/BIDdetect_data_invitro_delpos.txt` (gzip -9 -n) | 6abb0ec74cbba068248e47c875f82316 |
| `inputs/BIDdetect_data_incell_delpos.txt.gz` | `Figure3/plots/BIDdetect_data_incell_delpos.txt` (gzip -9 -n) | 9eacfaeb3b2dc4f33cc9eddf71f1e0db |
| `shipped/*.tsv` | `Figure3/plots/` paper outputs, unchanged | — |
| `legacy_rerun/invitro/` | `modification_analysis.R` (md5 8abf253f160547b5a046a9ef0e04aba3) rerun on the in-vitro input: `--prefix invitro_delpos --cores 1` | — |
| `legacy_rerun/incell/` | `incell_analysis.R` (md5 dd76753112f15a4746f1869b054df449) rerun on the in-cellulo input; only the user-configuration block changed (input/output paths, `num_cores <- 1`, `plot_all_sites <- FALSE`, `run_cell_type_analysis <- FALSE`) | — |

Reruns: Sherlock, R 4.3.2 with blme 1.0.6, lme4 1.1.35.5, car 3.1.3, dplyr 1.1.4, 2026-10-09
(`/scratch/users/rodell/NanoModAmp_build/wp4_legacy_ref/run_refs.sh`). The rerun is the
reference for code equivalence. The shipped tables come from slightly different script and package
versions (in vitro: same categories, p-values within 3e-4 relative); the test checks that their
categories / site sets are reproduced.

Findings while building these fixtures:

- **R-30 (PROPOSED):** both legacy scripts apply `p.adjust(..., "BH")` inside a `mutate()` on a
  tibble still grouped by `(chr, pos)`, so the "adjusted" p-values equal the raw ones (visible in
  the shipped tables: `p.value == p.adjust_diff`, `p.value == p.value.BH`). The port keeps this
  as `p_adjust = "legacy"` (default) and offers `p_adjust = "BH"`.
- readr writes subnormal doubles with a wrong exponent (e.g. 6.5187e-315 is written as
  `6.5187e-298`, PFKP_chr10_3112271 in `WT_mod_HepG2`); the test treats p < 1e-290 as equal.
- **Paper vs rerun drift (in cellulo WT_mod_Both):** rerunning the unchanged legacy script today
  reproduces the shipped categories at 759 of 760 sites; RPL22_chr1_6186768 is Inconclusive in the
  paper and Unmodified in the rerun (p 1.4e-5 vs 3.0e-22). 23 of 341 tested p-values differ by more
  than 1e-3 relative. The PUS7 union is identical (184 sites). The port equals the rerun exactly, so
  the drift comes from the original run's package versions / optimizer, not from the port.
