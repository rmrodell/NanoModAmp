# Published-pipeline harness: sources

Purpose: rerun the **published** scripts (PUS7regulation2026, commit `4fd285b`) on the golden package, so the new pipeline can be compared with both the legacy run (the scripts that actually produced the paper numbers, `dev/legacy/`) and the published code (Becca, 2026-10-06). Outputs: `golden/expected/published_rerun/` (reference only, not binding). Driver: `dev/golden/run_published.sh`.

`original/` holds byte-exact copies (`git show 4fd285b:<path>`). Working copies have only the edits listed in `ORIGINAL_DIFF.patch`: (a) keep the pre-dedup filtered BAM (`*_high_qual_primary.bam` → `retained/`, replacing `rm -rf "$TMP_DIR"`), (b) tool locations as overridable env vars (`MINIMAP2`, `SEQTK`, `UMICOLLAPSE`), `ml R` pinned to `R/4.2.0`, runtime `install.packages` replaced by an error, `analysis_endo.R` directory configurable plus the missing `/` (R-07), (c) `incell_analysis.R` user configuration read from env vars (`PUB_INPUT`, `PUB_OUTDIR`, `PUB_CORES`, `PUB_RUN_CELLTYPE`). Two R files have CRLF line endings in the paper repo; they are preserved.

| Working copy | Source (PUS7regulation2026 @ 4fd285b) | md5 (original, first 12) | Notes |
|---|---|---|---|
| endogenous/trim_map_dedup.sh | Figure2/file_prep/trim_map_dedup.sh | c63a1a354ecc | ONT linked trim (both orientations), umi_tools extract, `minimap2 -ax splice -uf`, `-F 2304`, `-q 30`, UMICollapse |
| mpra_invitro/trim_map_dedup_mpra.sh | Figure3/file_prep/trim_map_dedup_mpra.sh | 81e051035221 | as Figure 2 plus pool trim (`-m 120 -O 10`), `minimap2 -ax sr` |
| mpra_incell/trim_map_mpra_published.sh | **CONSTRUCTED** from Figure3/file_prep/trim_map_dedup_mpra.sh | — | see below |
| endogenous/trim_map_dedup_3prime_only.sh | **CONSTRUCTED** from Figure2/file_prep/trim_map_dedup.sh | — | 20250418 samples only; see below (D31, R-29) |
| common/BIDdetect.sh, common/bam_counts_fast.R | Figure2/BIDdetect/ (Figure 3 README: "the same scripts as in Figure2/BIDdetect") | fb5d8605c7c1, d5e1ae8f21fd | differs from the legacy counting scripts (kmer/ref-base logic in `bam_counts_fast.R`; empty-output guard in `BIDdetect.sh`) |
| common/sample_name.R | not in the paper repo; `rmrodell/BIDamplicon/legacy/sample_name.R` (D26), same file as `dev/legacy` | 48f9939d1482 | |
| mpra_invitro/modification_analysis.R | Figure3/mpra_sites/modification_analysis.R | 8abf253f1605 | byte-identical to the legacy copy (`4673e50`) |
| mpra_incell/incell_analysis.R | Figure3/mpra_sites/incell_analysis.R | dd76753112f1 | byte-identical to the legacy copy |
| endogenous/analysis_endo.R | Figure2/endo_sites/analysis_endo.R | 6145317f7424 | byte-identical to the legacy copy |

So the published and legacy runs share the site-calling code and differ in preprocessing and counting.

## Constructed in-cellulo variant (`mpra_incell/trim_map_mpra_published.sh`)

**Scope: this applies ONLY to the paper's golden in-cellulo samples (run 20241114).** That library was built without ONT adapters and without UMIs. The paper processed it with `dev/legacy/mpra_incell/trim_map_mpra.sh` (pool-adapter orientation passes, `minimap2 -a -k5`, no MAPQ filter, no UMI, no dedup). The published Figure 3 script requires ONT adapters and UMIs, so on these reads it would discard everything. **This combination never existed in the paper.** It is the closest "published" equivalent:

- orientation passes on the MPRA pool adapters, exactly as the legacy in-cellulo script: sense `GACGCTCTTCCGATCT...CACTCGGGCACCAAGGAC`, antisense `GTCCTTGGTGCCCGAGTG...AGATCGGAAGAGCGTC`, `-m 125 -O 15 --discard-untrimmed`, seqtk RC, concatenate;
- no `umi_tools extract`, no second pool trim (pool adapters are already removed), no UMICollapse;
- published mapping and filters kept: `minimap2 -ax sr`, `samtools view -F 2304`, `-q 30`; the final BAM (`final_bam/`) is the filtered primary MAPQ ≥ 30 BAM.

It is **not** a general in-cellulo protocol. The production pipeline's defaults keep ONT trimming + UMI extraction + dedup for all data, including in-cellulo data. The production pipeline reproduces this variant only when explicitly run with `orientation_adapters=pool, umi=false`, which the golden in-cellulo sample sheet/config sets for run 20241114.

## Constructed 3′-only endogenous variant (`endogenous/trim_map_dedup_3prime_only.sh`; D31, R-29)

**Scope: ONLY the golden endogenous run 20250418 samples.** Their reads lack the 5′ ONT adapter, so the published linked sense/antisense passes keep ~0 reads. The paper processed this run with a single 3′-adapter pass (`dev/legacy/endogenous/20250418/dedup_mapping.sh`). The variant replaces the published steps 2–3 (linked passes, RC, merge) with that pass: `cutadapt -m 125 -O 15 -a GAAGATAGAGCGACAGGCAAGT`, untrimmed reads kept (no `--discard-untrimmed`). Everything after is the unmodified published script: `umi_tools extract`, `minimap2 -ax splice -uf`, `-F 2304`, `-q 30`, UMICollapse. **This combination never existed in the paper.** Run 20251022 uses the unmodified published script. The production default still requires both adapters (D4); the 3′-only mode is an explicit opt-in (D31).

Test (293T_P36_1_BS, 650 package reads): 592 after 3′ trim, 494 primary MAPQ ≥ 30, 241 after UMICollapse (RHBDD2 115, HDAC6 126); the unmodified published script kept 0.

## Golden-run choices (run_published.sh)

- **Endogenous as one experiment.** `analysis_endo.R` reads one `BIDdetect_data.txt` and filters `vector == "WT"` / compares `WT` vs `KD`. Runs 20250418 and 20251022 are staged under normalised sample names: 20250418 `P36`, `pLKO` → `WT`, `P97`, `shPUS7` → `KD` (P97 = KD, Becca 2026-10-06); reps made run-unique (`0418r1`, `0418r3`, `1022r1`, `1022r2`) because rep numbers collide across runs. The mapping is saved as `name_map.tsv`. The package sample sheet for the new pipeline should use the same names.
- `analysis_endo.R` has no switch for its cell-type analyses; they run and are kept but are out of scope (Becca 2026-10-06). `incell_analysis.R` runs with `run_cell_type_analysis = FALSE`.
- In vitro: `modification_analysis.R --prefix invitro_delpos --plot_all_sites`, as in the paper invocation (Pool1/invitro/process.txt).
- Counting uses the selected-chr delpos BEDs (coordinates unchanged), as V1/V2.

## Findings from the test runs (2026-10-06)

- **Run 20250418 reads lack the 5′ ONT adapter.** In the test sample (650 package reads) the 5′ sense adapter occurs exactly in 2 reads and the 3′ sense adapter in 42; the published linked-adapter passes keep **0** reads. The paper's 20250418 processing (legacy `dedup_mapping.sh`) trimmed only the 3′ adapter. Resolved by D31 (Becca, 2026-10-06): the golden run uses the constructed 3′-only variant above for 20250418.
- UMICollapse is deterministic: two runs on the same input gave identical BAM records.
