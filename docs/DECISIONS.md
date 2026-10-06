# Decisions

Plan-level decisions (from `plan.md` §1; source Becca, 2026-09-28 unless dated):

| # | Topic | Decision |
|---|---|---|
| D1 | Users | Anyone, anywhere. Nothing specific to Sherlock or SLURM goes in pipeline code. |
| D2 | Entry point | Demultiplexed FASTQ (`.fastq.gz`, either one file or a directory of files per sample). |
| D3 | Plots | Standard pipeline output. Every plot is saved as **both PDF and PNG**. |
| D4 | Adapter trimming | Use cutadapt linked adapters, and keep **only reads with both adapters present**. This applies to every linked-adapter pass. |
| D5 | Mapping | minimap2 arguments are user-configurable. Defaults are set per library type: `mpra` = `-ax sr`, `endogenous` = `-ax splice -uf` (the paper settings). The docs must explain this choice and how to change it. |
| D6 | Endogenous reference | A set of transcript sequences, not a genome. |
| D7 | UMI dedup | Already validated. Keep it exactly as written (umi_tools extract with a 3′ 10-nt string pattern, then UMICollapse with default settings). |
| D8 | Length / overlap settings | The paper's MPRA-specific values become defaults, all configurable. |
| D9 | Denominator | Keep `totalReads = A+C+G+T+Deletion+Insertion` exactly as written, with a documentation note. Apply the minimum coverage filter at the earliest stage (counting). |
| D10 | Min coverage | Default is `totalReads > 20` (strict), as in the counting code. |
| D11 | BED | Either whole amplicons or single sites, whichever the user provides. |
| D12 | Bases reported | T-only by default. A QC mode reports all bases for background deletion analysis. |
| D13 | Model | `blme::bglmer` (the Figure 3 implementation) for **all** analyses, including endogenous. |
| D14 | Controls | Samples are explicitly labeled `treat = input` or `treat = BS` in the sample sheet. No implicit `noPUS → input` recoding. The golden sample sheet labels the paper's noPUS samples as `input`. |
| D15 | Thresholds | SESOI, the pre-filter, and FDR are kept exactly as written in the Figure 3 code. |
| D16 | Effect filter timing | Same as Figure 3: the Δ filter is applied **after** testing and BH. |
| D17 | Replicates | `rep` is paired by batch across conditions, so `(1|rep)` is meaningful. |
| D18 | Random effects | The random-effects structure stays exactly as written; extra terms are configurable per analysis. |
| D19 | Site key | `(chr, pos)`. Strand is always `+` (RNA amplicons). |
| D20 | Fidelity | Results must match the paper, except where a clearly documented Change Register entry explains the difference. |
| D21 | Golden data | Build a test-data package from Becca's local paper data and ship it **in the repo** as a QC test for new users. |
| D22 | Workflow manager | **Nextflow DSL2** using the nf-core template and conventions. |
| D23 | Environments | Per-process containers (Docker, Singularity/Apptainer), plus a conda profile as a fallback. |
| D24 | Analyses shipped | (a) Treatment analysis (BS vs input); (b) **generic factor tests** on any sample-sheet column and pair of levels, which generalize the PUS7-dependency and cell-type analyses. |
| D25 | Equivalence edge case | When mean input delrate = 0, the equivalence bound is infinite. **Keep as written** and document it as a known limitation. |
| D26 | `sample_name.R` | **Found** at `rmrodell/BIDamplicon/legacy/sample_name.R`. Use the original in the legacy harness only (WP5). The production pipeline uses the sample sheet instead. |
| D28 | Primer design | A standalone `nma-design` tool for endogenous amplicon and primer design, built as **WP9** in parallel with WP1–WP6. See `primer_plan.md`. |
| D27 | BED coordinate convention | **Unknown.** WP5 must work it out from the golden data (see R-15) and raise it at Gate G1. |
| D29 | Libraries without ONT adapters or UMIs (2026-10-06) | Two per-run parameters: `orientation_adapters` (default `ont`) and `umi` (default `true`). **The defaults always apply, for endogenous, MPRA in vitro and MPRA in cellulo data alike: ONT trimming, UMI extraction and deduplication.** `orientation_adapters=pool, umi=false` is an opt-in exception that exists only so the golden package can run the paper's in-cellulo samples (run 20241114), which were sequenced from a library built without ONT adapters or UMIs. It is set only in `conf/test_golden.config` for those samples, never by default or by `library_type` (R-26). |
| D31 | 3′-only ONT adapter libraries (2026-10-06) | Opt-in per-run parameter `ont_adapter_mode` (default `linked`: both ONT adapters required, D4). `three_prime_only` reproduces the legacy trim for libraries whose reads lack the 5′ ONT adapter: one pass `cutadapt -m {trim1_min_length} -O {trim1_min_overlap} -a <S3>`, untrimmed reads kept, no antisense pass or RC. It exists only for the golden endogenous samples from paper run 20250418 and is set only for them in `conf/test_golden.config`; the default applies to all other data (R-29). |
| D30 | Three-way golden comparison (2026-10-06) | The golden package is rerun with **both** the legacy scripts (what produced the paper numbers) and the published Figure 2/3 scripts. Both outputs are kept (`expected/legacy_rerun/`, `expected/published_rerun/`) and the new pipeline is compared against each in WP7. |

## Decision log

Every decision with its date and approver (plan §0 rule 5). Golden-package decisions G-D1–G-D8 are in `docs/plans/golden_test_package_plan.md` §1.

| Date | ID | Decision | Approved by |
|---|---|---|---|
| 2026-10-06 | G1-a/F1 | Legacy harness runs the scripts actually executed for the paper runs (run-dir copies matching run logs), not the published Fig2/3 `file_prep` scripts. Read IDs come from pre-dedup BAMs with no primary/MAPQ filter. | Becca |
| 2026-10-06 | G1-a/F2 | `endoG1G2.fasta` and its BEDs are used for all endogenous runs. Selected amplicon sequences verified identical to the per-run references. | Becca |
| 2026-10-06 | G1-a/F4 | Preflight finding F4 withdrawn: `pool1_chr.txt` names are valid pool1 oligo names; used as the MPRA chr list for both experiments. | Becca |
| 2026-10-06 | G1-b | P97 = KD (20250418 293T). Minimal sample sets approved as proposed: endogenous 20250418 reps 1,3 + 20251022 reps 1,2 (32); in vitro May run IV_noPUS_in_1, IV_noPUS_BS_2 (→ input), IV_PUS7_BS_1, IV_PUS7_BS_2 (4); in cellulo HepG2 reps 1,2 + 293T reps 1,3 × {P101,P102,P3,P4} × {input,BS} (32). | Becca |
| 2026-10-06 | G1-b | Cell-type-specific analyses are out of scope for the golden package; tests focus on treatment and PUS7 dependency. | Becca |
| 2026-10-06 | G1-b | Keep the 6 MPRA targets in `pool1_chr.txt` (no added oligos) and the 3 endogenous targets; uncovered categories (Inconclusive, input delrate = 0, Modified-not-PUS7-dependent) are accepted gaps. | Becca |
| 2026-10-06 | V2 env | Site-calling scripts run under R/4.3.2 with `~/R/x86_64-pc-linux-gnu-library/4.3`; `ggrepel` 0.9.8 installed there (all other packages already present). Counting stays on R/4.2.0. | Becca |
| 2026-10-06 | G1-d | Size policy replaces uniform `f` + 2% background: cap selected reads at **250 per (sample, chr)** (seeded, nested subsample; chrs below the cap kept whole) plus **150 background reads per sample** (seeded). Depth floor applies only to sites that start above it. If measured size > 25 MB, lower the cap (nested). V1 runs before the cap; V2 on the capped package is binding. | Becca |
| 2026-10-06 | G1-c | Harness edit (e): `umi_tools dedup --random-seed=${UMI_DEDUP_SEED:-20261004}` in the three dedup scripts (20250418, 20251022, in vitro). The paper runs set no seed, so their dedup representatives cannot be regenerated; seeding makes L1 reruns, V2 and Repro deterministic (verified: two seeded runs byte-identical). | Becca |
| 2026-10-06 | G1-c | V1 pass rule. **Exact** (all count columns, totalReads, delrate) for runs counted on non-deduplicated BAMs: endogenous 20250418 (pre-dedup `sorted/`) and in cellulo (no dedup). **Tolerance** for runs counted on dedup BAMs (endogenous 20251022, in vitro): per site × sample, `totalReads` exact (same dedup depth) and each `*.count` column within tolerance of the paper value. Provisional tolerance `|Δ| ≤ max(3, 5% of totalReads)`, `|Δdelrate| ≤ 0.03`; to be calibrated from the L1 paper-vs-rerun differences and reported. | Becca |
| 2026-10-06 | D29 | Opt-in per-run parameters `orientation_adapters` (default `ont`) and `umi` (default `true`). Defaults (ONT trim + UMI + dedup) apply to ALL data incl. in vitro and in cellulo. `pool`/`umi=false` is set only for the golden in-cellulo samples from paper run 20241114 (library built without ONT adapters/UMIs), in `conf/test_golden.config`. | Becca |
| 2026-10-06 | D30 | Golden package rerun with both the legacy scripts (`expected/legacy_rerun/`) and the published Figure 2/3 scripts (`expected/published_rerun/`); the new pipeline is compared against both in WP7. | Becca |
| 2026-10-06 | WP0 | Start WP0 on branch `wp/0-scaffold`, pushed to GitHub. | Becca |
| 2026-10-06 | D31 | Opt-in `ont_adapter_mode` (default `linked`, both ONT adapters required). `three_prime_only` (legacy 3′-only trim, untrimmed kept) is set only for the golden endogenous samples from run 20250418, whose reads lack the 5′ ONT adapter; the published rerun uses the same trim for them as a documented variant. | Becca |
| 2026-10-06 | WP0 | Nextflow minimum lowered from the template's 25.10.4 to 25.04.7, the newest Nextflow module on Sherlock; stub runs verified on 25.04.7. To confirm at G0. | Orchestrator (pending G0) |
