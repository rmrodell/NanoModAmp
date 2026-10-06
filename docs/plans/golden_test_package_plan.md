# Golden Test Package — Build Plan for Claude Agents

**Purpose:** Build a small, public test dataset from real paper data. It ships in `rmrodell/NanoModAmp` at `tests/data/golden/`. New users run it with `-profile test_golden` to confirm that their installation reproduces known results. CI uses the same dataset as a regression test.

This plan expands **WP5** of `plan.md`. Everything in `plan.md` §0 (agent rules) and §1 (decisions) still applies.

**Owner / reviewer:** Becca (rodell@stanford.edu). Steps marked **[GATE]** need her approval.

---

## 1. Decisions for this package (source: Becca, 2026-09-28)

| # | Topic | Decision |
|---|---|---|
| G-D1 | Experiments | Three: **endogenous** (Figure 2), **MPRA in vitro** (Figure 3), **MPRA in cellulo** (Figure 3) |
| G-D2 | Inputs provided | For each experiment, in separate directories: demultiplexed FASTQs, reference FASTA + BED, paper count tables, and paper final (deduplicated) BAMs |
| G-D3 | Read selection | By **chr values Becca provides** per experiment. Reads are assigned to a chr by re-running the **legacy preprocessing** and taking the read IDs whose **pre-dedup, primary, MAPQ ≥ 30** alignment is on a selected chr. The matching raw reads are then pulled from the original FASTQs, keeping all PCR duplicates. |
| G-D4 | Background reads | Add a seeded **~2% random sample** of all other raw reads per sample, so the adapter, length, mapping and MAPQ filters are exercised |
| G-D5 | Samples | The **minimal set** that can run every analysis for that experiment, with ≥ 2 batch-paired reps per condition |
| G-D6 | Size | Total package ≤ **25 MB**. If it's too big, **downsample reads** (seeded). Exact equality with the paper is then checked *before* downsampling (§5 V1), and the expected outputs come from a legacy rerun *after* downsampling (§5 V2). |
| G-D7 | Reference | Ship the **full reference** FASTA for each experiment, so mapping and MAPQ match the paper |
| G-D8 | Compute | Build on **Sherlock** using the original modules and scripts. The final package may be **public**. |

**Updates, 2026-10-06 (Becca; details in `docs/DECISIONS.md`).** These override the text below where they conflict:

- **Harness = scripts actually run (F1).** §4.1 copies the run-dir/commit versions that produced each paper run, not the published `file_prep` scripts; §4.3 read IDs come from the rerun pre-dedup BAM with no primary/MAPQ filter. Allowed edit (e): `umi_tools dedup --random-seed`.
- **Reference:** `endoG1G2.fasta` for all endogenous runs (F2). Selected targets: RHBDD2, HDAC6 (20250418), STIM1 (20251022); MPRA: the six `pool1_chr.txt` oligos.
- **Size (G1-d), replaces G-D4/G-D6 downsampling:** ≤ 250 selected reads per (sample, target) + 150 background reads per sample, seeded and nested.
- **V1 rule (G1-c):** exact for runs counted on non-deduplicated BAMs (20250418 pre-dedup, in cellulo); tolerance for dedup-counted runs (20251022, in vitro).
- **Two reruns (D30):** `expected/legacy_rerun/` (binding) and `expected/published_rerun/` (published Figure 2/3 scripts in `dev/published/`; the paper's in-cellulo samples via the constructed pool-adapter/no-UMI variant, D29). The pool/no-UMI settings apply only to these golden in-cellulo samples; the production pipeline's defaults keep ONT trimming, UMI extraction and dedup for all data.
- Cell-type-specific analyses are out of scope.
- **20250418 reads lack the 5′ ONT adapter (D31):** the published rerun and the new pipeline use the opt-in `ont_adapter_mode=three_prime_only` for those samples only; defaults require both adapters.

---

## 2. Inputs Becca provides

Becca fills in `dev/golden/golden_inputs.yaml`. A filled-in draft already exists at `STAR_Protocols/Code/golden_inputs.yaml`, and that file is the authoritative schema. Gate **[G1-a]** passes when no `TODO` values remain.

The layout is per run. Each experiment has experiment-level defaults and a `runs:` list:

```yaml
experiments:
  <experiment>:
    library_type: endogenous | mpra
    col_names: celltype_vector_rep_treat
    reference_fasta: ...            # default; a run may override it
    bed_amplicon: ...
    bed_position: ...               # optional
    bed_used_for_paper_counts: bed_amplicon | bed_position
    runs:
      - run_id: ...
        processing_script: fig2_trim_map_dedup | fig3_trim_map_dedup_mpra | legacy_dedup_mapping | legacy_dedup_pool | other:<path>
        fastq_dir: ...
        sample_map: ...
        paper_bam_dir: ...
        paper_counts: ...           # per-run table (endogenous) ...
    paper_counts: ...               # ... or one combined table for all runs (MPRA)
    chr_list: ...
    samples: auto
```

**Rules for multiple runs:**

- **Each run is preprocessed with its own `processing_script`, reference and sample map.** The legacy harness (§4.1) must support every script value used. `legacy_dedup_mapping` (BIDamplicon/legacy/dedup) differs from the Figure 2 script:
  - cutadapt trims only the 3′ adapter;
  - minimap2 runs with `-a -k5`;
  - there is no primary or MAPQ filter;
  - dedup uses umi_tools `--method directional`;
  - the deduplicated BAMs are in `minimap2/dedup/`, while `minimap2/sorted/` holds the pre-dedup BAMs.

  For those runs, check L1 against `minimap2/dedup/`, and take read IDs from the pre-dedup `sorted/` BAM, filtered the way that script filtered.
- **`sample_map: "none"`** (e.g., endogenous run 20250418) means the FASTQ directory already holds one file per sample. `sample_id` is then the file basename (`<sample>.fq`), as `dedup_map_folder.sh` used it.
  - Paper BAMs from `legacy_dedup_mapping` are named `<sample>_UMI_sort.bam` / `<sample>_UMI_dedup.bam`. `BIDdetect.sh` strips only `_sort`, so the paper `sample` field may end in `_UMI`. `sample_name.R` keeps only the first 4 fields, which drops that suffix.
  - The agent must reproduce this naming exactly when matching paper counts.
  - It must also determine whether the paper counts were made from the `sorted/` (pre-dedup) or the `dedup/` BAMs, by testing both against the counts. Report the answer at G1-c.
- **Reference and BED:** every endogenous run uses `endoG1G2.fasta` and the single-site BED (`bed_position`), per Becca.
- **`sample_id` must be unique across runs** within an experiment, or the plan must record how samples sequenced in several runs were merged. If any are duplicated, stop at G1-a.
- **Each run's counts are compared with its own paper table** when a table is given per run. Otherwise they are compared with the experiment's combined table.

**Preflight checks.** The agent runs these before any compute and fails with a clear message if any check fails:

0. There are no `TODO` values and no unexpanded `$` variables. Every run has `processing_script`, `fastq_dir`, `sample_map` and `paper_bam_dir`, and a `paper_counts` table at either the run or the experiment level. Some inputs stay on `/scratch` by choice; they won't be copied to `$OAK`. At preflight the agent records the md5 and modification time of every input in `provenance.json`. Before each later step, it checks that the scratch inputs still exist and are unchanged. If any is missing, it stops and asks Becca to restore it.
1. Every path exists and is readable.
2. Every chr in each `chr_list` is present in the reference FASTA index, the BED and the paper counts table. Report any that are missing.
3. Every sample in `sample_map` has FASTQ data.
4. Every paper BAM corresponds to a sample in `sample_map`.
5. The `col_names` split of the paper `sample` names reproduces the metadata columns in `paper_counts`.
6. The reference size is recorded. If the three references together exceed 15 MB compressed, **stop and ask**, because references can't be downsampled.

---

## 3. Selection

### 3.1 Sites (chr values)
These come from Becca's `chr_list` per experiment. The agent does **not** choose chrs. It reports a coverage profile of the list so Becca can confirm the list covers what the tests need. The profile shows, per chr, the paper categories of its sites (Modified, Unmodified, Inconclusive, PUS7-dependent, cell-type-specific), whether any site has input delrate = 0, and whether any site is near the coverage threshold. **[GATE G1-b]**

### 3.2 Samples (minimal set, G-D5)
When `samples: auto`, the agent proposes the smallest set of samples that lets every analysis in that experiment's section of `analyses_example.yaml` (from `plan.md` WP4) run with ≥ 2 reps per condition:

- **Endogenous:** WT and KD × input and BS × 2 reps, for each cell type used in the Figure 2 comparisons.
- **MPRA in vitro:** noPUS-BS (becomes `treat=input`, per D14) and PUS7-BS × 2 reps.
- **MPRA in cellulo:** HepG2 and 293T × {P102 (WT), P101 (KD), P4 (WT), P3 (OE)} × input and BS × 2 reps. This supports WT_mod, PUS7_dep and celltype_spec.

When the paper has more than 2 reps, pick reps 1 and 2 and record why. Becca approves the list **[GATE G1-b]**.

---

## 4. Build steps (Sherlock)

All scripts live in `dev/golden/`, and every step is idempotent (rerunning gives the same result). Submit through SLURM arrays that reuse the paper's `submit_file_prep.sh` pattern. Run every step inside `sherlock.workdir`, and never write to the input directories.

### 4.1 Legacy harness setup
- Copy the paper scripts (Figure 2 and Figure 3 `file_prep/`, `BIDdetect/`, `mpra_sites/`, `endo_sites/`) into `dev/legacy/`. Allowed edits:
  - **(a)** Keep the pre-dedup BAM. Replace `rm -rf "$TMP_DIR"` with a copy of `*_high_qual_primary.bam(.bai)` to `retained/`.
  - **(b)** Fix paths: tool locations, module names, and the missing `/` in `analysis_endo.R`.
  - **(c)** Fill in the `incell_analysis.R` user-configuration section.
  - **(d)** Add the original `sample_name.R` from `rmrodell/BIDamplicon/legacy/` (D26).
- Make **no logic edits.** Commit `dev/legacy/ORIGINAL_DIFF.patch`, which shows every change against the paper repo commit hash.
- Record the actual module versions loaded (`ml list`), plus `minimap2 --version`, the UMICollapse commit, `umi_tools --version` and `seqtk` in `dev/legacy/versions.txt`.

### 4.2 Legacy preprocessing on the selected samples (full reads)
- Run the legacy `trim_map_dedup*.sh` for each selected sample, using the paper's reference and parameters.
- **Check L1:** the read IDs in each rerun dedup BAM must equal the read IDs in the matching paper BAM, restricted to the selected chrs. Compare the sets, not the files byte for byte. If they don't match, **stop.** Report the size of the mismatch and likely causes (tool versions, non-determinism in UMICollapse) at **[GATE G1-c]**.

### 4.3 Read-ID assignment
- From each `retained/*_high_qual_primary.bam`, take the reads whose reference name is in `chr_list`. That BAM is already filtered to primary alignments with MAPQ ≥ 30.
- Normalize the IDs: strip the `_<UMI>` suffix that umi_tools appended, and strip any cutadapt/seqtk suffixes. Deduplicate the ID list.
- Write `selected_ids/<sample>.txt`, and record the counts per sample and chr.

### 4.4 Background reads
- From each sample's raw reads **not** in the selected IDs, take a seeded random sample of `background_fraction` (`seqkit sample -p 0.02 -s <seed>` on the ID complement).
- Write `background_ids/<sample>.txt`.
- These reads will map to other chrs or fail filters. That's intended: the shipped BED is restricted to the selected chrs, so the extra reads are processed but never counted.

### 4.5 Extract raw reads
- Concatenate the sample's raw FASTQ(s) in the same order the legacy script used (lexical `find` order; record it).
- Use `seqkit grep -f` with the union of selected and background IDs to extract raw records **in their original orientation with adapters intact**.
- Write `stage/<experiment>/fastq/<sample_id>.fastq.gz` (a single file per sample).
- **Check L2:** the number of extracted reads equals the number of IDs. Any ID not found in the raw FASTQ is a hard error.

### 4.6 Lossless check, then size and downsampling
1. **Check V1 (lossless selection).** Run the legacy pipeline (preprocess → count → sample_name) on `stage/` with the **full reference** and the **selected-chr BED** (§4.7). For every selected site and sample, the counts must equal the paper counts table **exactly** (`totalReads`, all `*.count` columns, `delrate`). This proves the selection loses nothing. If it fails, stop at **[GATE G1-c]**.

   **Why an exact match is expected here, even though the package is smaller than the paper dataset.** V1 compares *per-sample, per-site counts*, and it runs *before* downsampling. Those counts depend only on the reads of that sample that reach that site:

   - **Choosing a subset of samples or replicates (G-D5) doesn't matter.** Each sample is processed independently through preprocessing, dedup and counting. Dropping replicate 3 cannot change the counts for replicate 1.
   - **Keeping only the selected chrs doesn't matter.** Every read that the paper run placed on a selected chr is included (§4.3). Dedup groups reads by alignment position and UMI, so every duplicate group at those positions is complete. Reads on other chrs never enter those groups.
   - **The same reads go through the same steps.** Each read's fate through trimming, UMI extraction, mapping against the full reference and MAPQ filtering is decided per read, so it is the same in the subset as in the full run. `seqkit grep` keeps the original read order, so tie-breaking in dedup is unchanged too.
   - **Background reads (§4.4) don't matter.** By construction, none of them ended on a selected chr in the full run, and the pipeline is deterministic, so they can't do so in the rerun.

   V1 covers counts only. It does **not** cover site calls. Site calls *are* expected to differ from the paper, for three reasons. Fewer samples and replicates change the model fits. Fewer sites change the Benjamini–Hochberg correction set. After downsampling, lower depth changes both the counts and the calls. That is why the binding expected outputs come from V2, not from the paper (§4.8).

   If V1 fails, the likely cause is not the selection but a difference between the rerun and the paper run. Examples: different tool versions; paper BAMs built from a different or combined set of runs (e.g., `combineruns.sh` in BIDamplicon/legacy); or FASTQs that don't match the paper BAMs. Check L1 (§4.2) should catch these first.
2. **Measure the size.** Add up all FASTQs (gz), references (gz), BED, tables and expected outputs, for each experiment and in total.
3. **If the total is ≤ `size_cap_mb`**, skip the rest of this step and go to 4.7.
4. **Otherwise, downsample:**
   - Find a single fraction `f` that brings the total under the cap. Apply it with seed `seed` to each sample's **selected** reads, stratified by chr, so every chr is reduced by the same fraction. Use `seqkit sample -p f -s seed` per (sample, chr) ID list.
   - Background reads are downsampled by the same `f`.
   - **Depth floor:** after downsampling, estimate the dedup depth per selected site (paper dedup depth × f). If any site falls below `min_dedup_depth_after_downsample`, first try a per-experiment `f` instead of a global one. If the floor still can't be met under the cap, **stop and report** the sizes per chr, and let Becca drop chrs or raise the cap **[GATE G1-d]**.
   - Record `f`, the seed and the per-sample, per-chr read counts before and after in `provenance.json`.
   - Because downsampling breaks exact equality with the paper, the paper numbers become a *reference* (§5). The binding expected outputs come from V2.

### 4.7 Reference, BED, sample sheet and analyses
For each experiment:

- `reference/<name>.fa.gz` (full, per G-D7) plus `.fai`.
- `targets.bed`: the paper BED, filtered to the selected chrs, with **coordinates unchanged**.
- `samplesheet.csv` in the `plan.md` §5.1 format:
  - `fastq` paths are relative to the package.
  - `treat` is explicit. MPRA in vitro noPUS-BS samples get `treat=input` (D14), noted in a `notes` column.
  - Metadata columns come from `col_names`.
  - `rep` is batch-paired.
- `analyses.yaml`: the experiment's part of `assets/analyses_example.yaml`, restricted to the included samples.

### 4.8 Expected outputs
Write these under `expected/` for each experiment:

- **`paper_reference/`**: rows from the paper counts table, and the paper site-calling tables, for the selected chrs and samples. For MPRA, use the shipped paper tables (`invitro_modification_significance.tsv`, `WT_mod_Both_significance.tsv`, `PUS7_dep_union_significant_summary.tsv`, etc.). For endogenous, use Becca's Figure 2 outputs if she has them. These are **reference only, never binding**:
  - The paper **counts** should equal `legacy_full_selected` (V1). They will not equal the downsampled package.
  - The paper **site calls** came from all samples, all replicates and all sites. They are expected to differ from calls on the package, because of different model inputs, a different BH set and lower depth. The WP7 concordance report shows the agreement for information only.
- **`legacy_full_selected/`**: the V1 legacy outputs (counts) on the non-downsampled selection. Include this only if it fits the cap; otherwise store md5s and summary statistics only.
- **`legacy_rerun/`** (V2, **binding**): run the full legacy chain on the **final** package, i.e. preprocessing, `BIDdetect.sh`, and the legacy site calling (`modification_analysis.R` for in vitro, `incell_analysis.R` for in cellulo, `analysis_endo.R` for endogenous). Store:
  - The read funnel per sample.
  - `BIDdetect_data.txt`.
  - All site-calling tables.
  - md5s of the dedup BAMs (not the BAMs themselves).
  - The plots, as PNG only, to save space.
- **`published_rerun/`** (D30, reference): the same chain with the published Figure 2/3 scripts (`dev/published/`), stored in the same layout as `legacy_rerun/`.
- Endogenous legacy calls use `glmer` (Figure 2 method), while the new pipeline uses `bglmer` (D13, R-21). So for endogenous, `legacy_rerun/` calls are **reference only**, and counts are binding. State this in the README.

### 4.9 BED convention finding (R-15)
Use the V1 data to find the convention. For each BED region, compare `start` with the smallest `pos` reported by the paper counts, and look for positions at `start` (0-based) that shouldn't exist in single-site regions. Leading hypothesis: the endogenous BEDs came from the BIDamplicon design tool, which writes **standard 0-based** BEDs (see `primer_plan.md` §6). If so, the missing `+1` in `bam_counts_fast.R` is a bug. Write `dev/golden/R15_bed_convention.md` with the evidence (tables of examples) and a recommendation for `bed_coordinates`. **[GATE G1-e]**

### 4.10 Package and manifest
The final layout:

```
tests/data/golden/
  README.md                 # provenance, decisions, how to run, what's binding vs reference
  provenance.json           # paper repo commit, input paths (redacted to basenames), tool versions,
                            # seed, f, per-sample/per-chr read counts (selected, background, before/after)
  MANIFEST.md5              # md5 of every file in the package
  endogenous/ mpra_invitro/ mpra_incell/
    fastq/<sample_id>.fastq.gz
    reference/<ref>.fa.gz(.fai)
    targets.bed
    samplesheet.csv
    analyses.yaml
    expected/{paper_reference, legacy_full_selected (if kept), legacy_rerun}/
```

- Compress FASTQs with `gzip -9`, or with bgzip if that's smaller. Record the choice.
- Before committing, verify the total size is ≤ `size_cap_mb`.
- Also add `conf/test_golden.config`, pointing at the three sample sheets. Running all three may need three runs or a param. Coordinate this with the Orchestrator and record it in the `plan.md` contracts changelog.

---

## 5. Validation summary

| Check | What | Pass criterion | On failure |
|---|---|---|---|
| Preflight | Paths, chr lists, sample maps, name splitting | All pass | Stop, report |
| L1 | Rerun dedup BAM read IDs vs paper BAM, selected chrs | Identical sets | **[GATE G1-c]** |
| L2 | Extracted reads = ID list | Exact | Hard error |
| V1 | Legacy **counts** on the selection, *before downsampling*, vs paper counts (site calls not compared) | Exact, every selected site × included sample | **[GATE G1-c]** |
| Size | Package total | ≤ cap after downsampling, depth floor met | **[GATE G1-d]** |
| V2 | Legacy rerun on the final package | Completes; becomes the binding expected output | Fix the harness |
| Repro | Rebuild from `golden_inputs.yaml` with the same seed | Byte-identical `MANIFEST.md5` | Fix non-determinism |
| Smoke | `nextflow run . -profile test_golden,docker` once the new pipeline exists | Runs; concordance report generated (`plan.md` WP7) | Handed to WP7 |

---

## 6. Deliverables and definition of done

- `dev/golden/` contains the build scripts (`build_golden.sh` plus Python/R helpers), `golden_inputs.template.yaml`, SLURM submit wrappers and `R15_bed_convention.md`.
- `dev/legacy/` contains the patched scripts, `ORIGINAL_DIFF.patch`, `versions.txt` and the reconstructed `sample_name.R`.
- `tests/data/golden/` contains the complete package (§4.10), within the size cap.
- `PROGRESS.md` is updated, and every gate outcome is recorded in `docs/DECISIONS.md`.
- **Done when** L1, L2, V1, V2 and Repro all pass, Becca has approved G1-a through G1-e, and the package is committed via PR to `dev`.

## 7. Gates

| Gate | Becca provides / decides |
|---|---|
| G1-a | Completed `golden_inputs.yaml` (directories, sample maps, `col_names`, chr lists) |
| G1-b | Approves the chr coverage profile and the proposed minimal sample list |
| G1-c | Only if L1 or V1 fails: decides how to proceed (e.g., accept documented tool-version drift) |
| G1-d | Only if the depth floor can't be met under the cap: drops chrs or raises the cap |
| G1-e | Picks the BED coordinate convention (R-15) based on the evidence |
