# NanoModAmp — Execution Plan for Claude Agents

**Goal:** Build a production pipeline for Nano-BID-Amp data. It goes from demultiplexed FASTQ to site calls and standard plots, and it must run anywhere: a laptop, any HPC scheduler, or the cloud. The pipeline must reproduce the results of *RNA sequence, structure, and cell type specific features drive pseudouridylation by PUS7* (Figures 2 and 3). The only allowed differences are ones explained by a documented entry in the Change Register (§9).

- **Target repo (write here):** `https://github.com/rmrodell/NanoModAmp`. It is currently a skeleton, and all its files are empty. Put a copy of this `plan.md` at the repo root.
- **Reference repo (read-only):** `https://github.com/martinezlab/PUS7regulation2026`. Use only `Figure2/` and `Figure3/`.
- **Owner / reviewer:** Becca (rodell@stanford.edu). Any decision marked **[GATE]** needs her approval.

---

## 0. Rules for every agent

1. **Scientific behavior is frozen.** Port the algorithms in §6 exactly. Change behavior only for a Change Register entry (§9) that has status `APPROVED`. If you find a new discrepancy, add it to the register as `PROPOSED` and raise it at the next gate. Do not fix it on your own.
2. **Contracts are frozen after WP0.** The file formats in §5 are the interfaces between parallel work packages. Only the Orchestrator can change them, and every change goes in `docs/contracts/CHANGELOG.md`.
3. **No hidden environment.** No hard-coded paths, no `module load`, no installing packages at runtime, no `detectCores()` defaults. Every tool runs in a pinned container. Every parameter lives in `nextflow_schema.json` or the analysis config.
4. **Tests before merge.** A work package is done only when its acceptance criteria (§7) pass in CI.
5. **Track progress.** Update `PROGRESS.md` (WP status, blockers, open questions) whenever you finish or get blocked. Log every decision in `docs/DECISIONS.md` with its date and who approved it.
6. **Git workflow.** Work on `wp/<id>-<short-name>` branches, one PR per work package, with the Orchestrator merging into `dev`. `main` receives only gate-approved releases.
7. **When unsure, stop.** Add the question to `PROGRESS.md` under "Open questions" and continue with unblocked work. Don't guess on scientific choices.

---

## 1. Decisions already made (source: Becca, 2026-09-28)

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
| D27 | BED coordinate convention | **Decided at G1-e (2026-10-06):** default `bed_coordinates=bed0` (standard 0-based, half-open). `one_based_start` reproduces the legacy reading (start and end both 1-based, inclusive) for the paper's BEDs. Under `bed0`, a row with start = end is an error pointing to `one_based_start` (no auto-detect). The golden `targets.bed` is converted to standard 0-based (start = site − 1, end = site). See `dev/golden/R15_bed_convention.md`. |
| D29 | Libraries without ONT adapters or UMIs (2026-10-06) | Two per-run parameters: `orientation_adapters` (default `ont`) and `umi` (default `true`). **The defaults always apply, for endogenous, MPRA in vitro and MPRA in cellulo data alike: ONT trimming, UMI extraction and deduplication.** `orientation_adapters=pool, umi=false` is an opt-in exception that exists only so the golden package can run the paper's in-cellulo samples (run 20241114), which were sequenced from a library built without ONT adapters or UMIs. It is set only in `conf/test_golden.config` for those samples, never by default or by `library_type` (R-26). |
| D31 | 3′-only ONT adapter libraries (2026-10-06) | Opt-in per-run parameter `ont_adapter_mode` (default `linked`: both ONT adapters required, D4). `three_prime_only` reproduces the legacy trim for libraries whose reads lack the 5′ ONT adapter: one pass `cutadapt -m {trim1_min_length} -O {trim1_min_overlap} -a <S3>`, untrimmed reads kept, no antisense pass or RC. It exists only for the golden endogenous samples from paper run 20250418 and is set only for them in `conf/test_golden.config`; the default applies to all other data (R-29). |
| D32 | Merging runs at the counts level (2026-10-06) | Runs that need different preprocessing parameters (e.g. golden endogenous 20250418 with `ont_adapter_mode=three_prime_only` and 20251022 with defaults) are preprocessed and counted separately, then merged at the **counts** level for site calling. New parameter `input_counts`: one or more `counts_merged.tsv` files (§5.3) from earlier runs; when set, preprocessing and counting are skipped, the tables are merged by column name (R-06), and calling runs on the merged table. Used by `test_golden` for endogenous. |
| D33 | Endogenous random effects (2026-10-09) | Replicate and cell type: `(1|rep)` is always in the model (§6.3/§6.4); endogenous analyses that combine both cell types add `(1|celltype)`. Single-cell-type analyses use `(1|rep)` only. |
| D30 | Three-way golden comparison (2026-10-06) | The golden package is rerun with **both** the legacy scripts (what produced the paper numbers) and the published Figure 2/3 scripts. Both outputs are kept (`expected/legacy_rerun/`, `expected/published_rerun/`) and the new pipeline is compared against each in WP7. |

---

## 2. Architecture

```
samplesheet.csv ─┐
params / config ─┤
                 ▼
 ┌──────────────── PREPROCESS (per sample) ────────────────┐
 │ CAT_FASTQ → CUTADAPT_SENSE ┐                            │
 │            → CUTADAPT_ANTISENSE → SEQTK_RC ┘→ MERGE     │
 │ → UMITOOLS_EXTRACT → [mpra: CUTADAPT_POOL]             │
 │ → MINIMAP2 → SAMTOOLS_SORT → FILTER(-F 2304, -q 30)    │
 │ → UMICOLLAPSE → SAMTOOLS_INDEX → READ_FUNNEL metrics   │
 └────────────────────────────────────────────────────────┘
                 ▼
 COUNT (per sample, R/Rsamtools) → per-sample counts.tsv
                 ▼
 MERGE_COUNTS (+ sample-sheet metadata, join by name) → counts_merged.tsv
                 ▼                          └→ [qc mode] BACKGROUND_QC
 CALL (per analysis in analyses.yaml, R package)
   treatment | factor  → tables + PDF/PNG plots
                 ▼
 SITE_SETS (optional union/intersection)  →  MULTIQC + run report
```

**Languages:**

- Nextflow handles orchestration.
- The R package `nanomodamp` (in `R/nanomodamp/`, installed in the `nanomodamp-r` container) handles counting, merging, statistics and plots. It is called through thin CLIs in `bin/`: `nma_count.R`, `nma_merge.R`, `nma_call.R`, `nma_sitesets.R`, `nma_bgqc.R`.
- Python (`tests/simulate/`, `dev/`) is used for the simulator and the test-package builder only.

**Repo layout (target):**

```
main.nf  nextflow.config  nextflow_schema.json  modules.json
workflows/nanomodamp.nf
subworkflows/local/{preprocess,count,call}.nf
modules/nf-core/...           # cutadapt, seqtk/seq, umitools/extract, minimap2/align,
                              # samtools/{sort,view,index}, umicollapse, multiqc
modules/local/{count,merge_counts,call_sites,site_sets,background_qc,read_funnel}/
conf/{base,modules,test,test_golden,test_full}.config
assets/{schema_input.json, schema_analyses.json, adapters_default.yaml,
        analyses_example.yaml, multiqc_config.yml}
bin/nma_*.R
R/nanomodamp/                 # DESCRIPTION, R/, tests/testthat/
containers/nanomodamp-r/Dockerfile   (+ legacy/ for WP5)
tests/
  simulate/                   # simulator (Python) + truth generation
  data/synthetic/             # small generated dataset + truth tables (committed)
  data/golden/                # test package from paper data (committed, ≤ 25 MB)
  nf-test/                    # nf-test specs for modules, subworkflows, pipeline
dev/                          # golden/ (package builder), legacy/ (scripts as run for the paper),
                              # published/ (published Fig 2/3 scripts); not shipped to users
docs/{usage.md, output.md, parameters.md, methods.md, contracts/, DECISIONS.md}
CHANGE_REGISTER.md  PROGRESS.md  CHANGELOG.md  CITATION.cff
.github/workflows/{ci.yml, linting.yml}
```

---

## 3. Parameters (defaults are the paper values)

| Param | Default | Notes |
|---|---|---|
| `input` | — | Sample sheet CSV (§5.1) |
| `library_type` | — | `endogenous` or `mpra`, one per run |
| `fasta` | — | Transcripts (endogenous) or oligo pool (MPRA). `.fai` is built if missing. |
| `bed` | — | ≥ 6 columns; amplicons or single sites (D11); strand should be `+` (warn otherwise) |
| `bed_coordinates` | `bed0` | `bed0` (standard, default) or `one_based_start` (legacy reading of the paper's BEDs). Under `bed0`, start = end rows are an error. See R-15, D27 |
| `ont_adapter_sense_5p` | `TTTCTGTTGGTGCTGATATTGCG` | |
| `ont_adapter_sense_3p` | `GAAGATAGAGCGACAGGCAAGT` | |
| `ont_adapter_antisense_5p` | `ACTTGCCTGTCGCTCTATCTTC` | Reverse complement of sense 3′ — validate at startup |
| `ont_adapter_antisense_3p` | `CGCAATATCAGCACCAACAGAAA` | Reverse complement of sense 5′ — validate at startup |
| `trim1_min_length` | `125` | cutadapt `-m` |
| `trim1_min_overlap` | `15` | cutadapt `-O` |
| `pool_adapter_5p` | `GACGCTCTTCCGATCT` | MPRA only |
| `pool_adapter_3p` | `CACTCGGGCACCAAGGAC` | MPRA only |
| `trim2_min_length` | `120` | MPRA only |
| `trim2_min_overlap` | `10` | MPRA only |
| `require_both_adapters` | `true` | D4. `false` reproduces legacy behavior (for diagnostics only) |
| `ont_adapter_mode` | `linked` | D31. `linked`: sense/antisense passes require both ONT adapters (D4). `three_prime_only` (opt-in, golden 20250418 samples only): single `-a <S3>` pass, untrimmed reads kept, no antisense pass/RC. Only valid with `orientation_adapters=ont` |
| `orientation_adapters` | `ont` | D29. `ont`: orientation passes use the ONT adapters (then the MPRA pool trim, if `mpra`). `pool`: orientation passes use the MPRA pool adapters with `trim1_*` settings, and the separate pool trim is skipped. `pool` requires `library_type=mpra`. |
| `pool_adapter_antisense_5p` | `GTCCTTGGTGCCCGAGTG` | Used only with `orientation_adapters=pool`. Reverse complement of `pool_adapter_3p` — validate at startup |
| `pool_adapter_antisense_3p` | `AGATCGGAAGAGCGTC` | Used only with `orientation_adapters=pool`. Reverse complement of `pool_adapter_5p` — validate at startup |
| `umi` | `true` | D29. `false` skips UMI extraction and UMICollapse; the filtered BAM is the final BAM. Funnel steps `umi_extracted` and `dedup` are omitted. |
| `umi_pattern` | `NNNNNNNNNN` | umi_tools `--extract-method=string --3prime`. Ignored when `umi=false` |
| `minimap2_args` | by library type | `mpra: -ax sr`; `endogenous: -ax splice -uf` (D5) |
| `min_mapq` | `30` | `samtools view -q` (≥ 30) |
| `sam_exclude_flags` | `2304` | Secondary + supplementary |
| `umicollapse_args` | `''` | Defaults (D7) |
| `min_coverage` | `20` | Keep sites with `totalReads > min_coverage` (D10), applied at counting |
| `pileup_max_depth` | `200000` | Warn if any position reaches it |
| `pileup_min_mapq` / `pileup_min_base_quality` | `1` / `1` | As in paper |
| `count_all_bases` | `false` | D12 QC mode: emits all-base counts + background report |
| `kmer_size` | `5` | Centered; NA near ends |
| `analyses` | — | YAML file of analyses (§5.4); if absent, counting only |
| `input_counts` | — | D32. Comma-separated `counts_merged.tsv` paths (or a directory of them). Skips preprocessing/counting; merges by name, then calls. Mutually exclusive with `input`. Errors: duplicate `sample_id` across tables; differing count/site columns. Metadata columns are the union (missing → NA, with a warning). Sites missing in a sample are simply absent (no fill). |
| `plot_all_sites` | `false` | Per-analysis override allowed |
| `allow_region_failures` | `false` | If false, any per-region counting error fails the task (R-12) |
| `outdir` | `results` | |

Profiles: `docker`, `singularity`, `apptainer`, `conda`, `test` (synthetic), `test_golden` (paper subset), `test_full` (documented, not in CI). Executor settings (SLURM, AWS Batch, etc.) are handled by the user's Nextflow config. The docs include an example SLURM config and one for Sherlock.

---

## 4. Output layout

```
results/
  pipeline_info/            # execution report, timeline, trace, software_versions.yml, params.json
  preprocess/<sample>/      # tool logs (cutadapt, umi_tools, minimap2, umicollapse)
  bam/<sample>.bam(.bai)    # final deduplicated BAMs
  metrics/read_funnel.tsv   # all samples, all steps (§5.2)
  counts/per_sample/<sample>.counts.tsv
  counts/counts_merged.tsv  # BIDdetect_data.txt equivalent (§5.3)
  qc/background/            # only if count_all_bases (tables + PDF/PNG)
  calling/<analysis_name>/
      <name>_log.txt  <name>_resolved_config.yaml
      data_summary/   data_raw/   plots/   # plots in PDF + PNG
  site_sets/<set_name>/     # optional unions/intersections
  multiqc/multiqc_report.html
```

---

## 5. Contracts (frozen after WP0)

### 5.1 Sample sheet (`assets/schema_input.json`, validated with nf-schema)
Required columns:

- `sample_id`: unique; matches `[A-Za-z0-9._-]+`.
- `fastq`: a `.fastq.gz` file, **or** a directory whose `*.fastq.gz` files are concatenated in lexical order.
- `treat`: exactly `input` or `BS`.
- `rep`: integer or string, the batch-paired replicate (D17).

Any other columns (e.g., `celltype`, `vector`, `batch`) are **metadata**. They are carried into the merged counts and can be used for analysis subsets and factors.

Validation errors:

- Duplicate `sample_id`.
- Missing file.
- An empty FASTQ directory.
- A `treat` value other than `input`/`BS`.
- Blank lines are ignored with a warning.

### 5.2 Read funnel (`metrics/read_funnel.tsv`)
Columns: `sample_id, step, file, records, bytes`.

Steps, in order:

1. `raw` (count `.gz` files correctly; R-08)
2. `trim1_sense`, `trim1_antisense`, `antisense_rc`, `merged`
3. `umi_extracted`
4. `trim2_pool` (MPRA only)
5. `mapped`, `primary`, `mapq_filtered`
6. `dedup`

Steps that a run's parameters skip (`umi_extracted` and `dedup` when `umi=false`; `trim2_pool` when `orientation_adapters=pool`) are omitted, not written as zero.

Also record `n_reads_in_both_orientations` (read IDs found in both the sense and antisense outputs).

### 5.3 Counts
Per-sample file `<sample>.counts.tsv`. Column names and order are **exactly** as follows, to stay compatible with the paper:

```
chr pos gene totalReads A.count C.count G.count T.count Deletion.count Insertion.count ref kmer strand delrate
```

`counts_merged.tsv` = `sample_id` + all sample-sheet metadata columns (in sample-sheet order) + the per-sample columns above. With `input_counts` (D32), several such tables are merged into one with the same layout (metadata columns in first-seen order), written as `counts/counts_merged.tsv` plus `counts/merge_sources.tsv` (sample_id → source table). Tables are **joined by column name, never by position** (R-06). Types: `pos` is an integer, and the count columns are integers.

### 5.4 Analyses config (`assets/schema_analyses.json`)
```yaml
analyses:
  - name: incellulo                 # treatment analysis
    type: treatment
    subset: {celltype: [HepG2, 293T]}      # optional; column -> allowed values
    random_effects: ""            # extra terms appended to "delrate ~ treat + (1|rep)"
    sesoi: 0.05
    fdr: 0.05
    plot_all_sites: true
    colors: {modified: "#c154c1", input: "#eee8aa"}
  - name: PUS7_dep_Both_WT_v_KD   # generic factor analysis
    type: factor
    subset: {celltype: [HepG2, 293T], vector: [P101, P102]}
    factor: vector
    levels: [P101, P102]          # [baseline, experimental]; dd = experimental - baseline
    random_effects: "(1|celltype)"
    direction: positive           # positive | both
    sesoi: 0.05
    fdr: 0.05
  - name: PUS7_dep_Both_WT_v_OE
    type: factor
    subset: {celltype: [HepG2, 293T], vector: [P4, P3]}
    factor: vector
    levels: [P4, P3]              # P4 = WT (baseline), P3 = OE (experimental)
    random_effects: "(1|celltype)"
site_sets:                        # optional
  - name: PUS7_dep_union
    op: union                     # union | intersection | difference
    of: [PUS7_dep_Both_WT_v_KD, PUS7_dep_Both_WT_v_OE]
    quartiles_by: {analysis: incellulo, column: delta_delrate}   # optional
```
`assets/analyses_example.yaml` must contain the full set of analyses that reproduces the paper's Figure 3 in-cellulo, Figure 3 in-vitro and Figure 2 results (built in WP4/WP5).

### 5.5 Site-calling result tables
- **Treatment**, `data_summary/<name>_modification_significance.tsv`. Columns: `chr, pos, category, avg_delrate_BS, avg_delrate_input, delta_delrate, p.adjust_diff, equivalence_status, site_specific_eqbound`, sorted by `delta_delrate` descending. Additional columns may be appended: `p.value`, `is_equivalent`, `model_status`.
- **Factor**, `data_summary/<name>_significant_summary.tsv`. Columns: `chr, pos, term, p.value, note, p.value.BH, delta_delrate_<baseline>, delta_delrate_<experimental>, dd_delrate`. Also write the all-sites table `<name>_all_tests.tsv`, which includes non-significant sites and failed sites with their `note`.

---

## 6. Algorithm specifications (port exactly)

### 6.1 Preprocessing (source: `Figure3/file_prep/trim_map_dedup_mpra.sh`, `Figure2/file_prep/trim_map_dedup.sh`)
Steps 2–5 and 8 depend on `orientation_adapters` and `umi` (D29):

| Mode | Orientation passes (step 2) | Step 4 UMI | Step 5 pool trim | Step 8 dedup |
|---|---|---|---|---|
| `ont` + `umi=true` (**default for all data**, including MPRA in cellulo) | ONT adapters | yes | MPRA only | UMICollapse |
| `pool` + `umi=false` (opt-in; only the golden in-cellulo samples from paper run 20241114) | pool adapters `<P5>...<P3>` and `<PA5>...<PA3>`, `trim1_*` settings | skipped | skipped (already trimmed) | skipped |
| `ont` + `umi=false` | ONT adapters | skipped | MPRA only | skipped |

`pool` + `umi=true` is rejected at startup (the UMI sits outside the pool adapters and would be trimmed away).

With `ont_adapter_mode=three_prime_only` (D31, opt-in), step 2 is a single `cutadapt -m {trim1_min_length} -O {trim1_min_overlap} -a <S3>` pass that keeps untrimmed reads; step 3 (RC + merge) is skipped and the funnel records `trim1_3prime` instead of the sense/antisense steps.

1. **Input:** a single file, or a lexically sorted concatenation of the directory's `*.fastq.gz` files.
2. **Sense pass:** `cutadapt -m {trim1_min_length} -O {trim1_min_overlap} --discard-untrimmed -g "<S5>;required...<S3>;required"`. Run the antisense pass the same way with `<A5>...<A3>`. Use the `;required` modifiers, or whatever the pinned cutadapt version needs, so that **both** adapters are required (D4, R-04). WP2 must prove this with a unit test.
3. Reverse-complement the antisense reads (`seqtk seq -r`), then concatenate the sense reads followed by the RC antisense reads. Keep reads that appear in both outputs, as the paper did, but count them in the funnel. If they exceed 0.1% of reads, flag it as an open question.
4. `umi_tools extract --extract-method=string --bc-pattern={umi_pattern} --3prime`. The UMI goes into the read name with the `_` separator that UMICollapse expects.
5. **MPRA only:** `cutadapt --discard-untrimmed -m {trim2_min_length} -O {trim2_min_overlap} -g "<P5>;required...<P3>;required"`.
6. `minimap2 {minimap2_args} -t {task.cpus} {fasta} reads.fq`, then `samtools sort`.
7. `samtools view -b -F {sam_exclude_flags}`, then `samtools view -b -q {min_mapq}`. The two steps can be merged into one call because the result is identical. Then index.
8. `umicollapse bam -i in.bam -o out.bam {umicollapse_args}`, then index.
9. The md5 copy step from the legacy scripts is replaced by Nextflow's `publishDir`. No md5 step is needed.

### 6.2 Counting (source: `Figure2/BIDdetect/bam_counts_fast.R`)
For each BED region (chr, start, end, gene, strand):

1. **Region start:** apply the convention chosen by `bed_coordinates`. Legacy used the BED start as-is (see R-15).
2. **Pileup:** `Rsamtools::pileup` with `PileupParam(max_depth, min_mapq=1, distinguish_nucleotides=TRUE, ignore_query_Ns=TRUE, min_base_quality=1, include_insertions=TRUE, include_deletions=TRUE, distinguish_strands=TRUE)`.
3. Keep only pileup rows on the region's strand. For `-` strand, complement the nucleotides with `chartr("ATCG","TAGC")`.
4. Reshape wide, summing counts by `(seqnames, pos, nucleotide)`. `-` becomes `Deletion` and `+` becomes `Insertion`.
5. `ref` = the reference base, complemented on `-` strand. `kmer` = the 5-mer at `pos−2..pos+2`, reverse-complemented on `-` strand. `kmer` is NA if `pos < 3` or `pos > chrom_len − 2`.
6. Unless `count_all_bases` is set, keep only `ref == "T"`.
7. Fill missing base columns with 0. Set `totalReads = A+C+G+T+Deletion+Insertion` (D9; note that insertions are double-counted). Keep `totalReads > min_coverage`. Set `delrate = Deletion/totalReads`.
8. Sites that appear in several overlapping BED regions: output the rows as the legacy code did, and record the count of duplicates in the log.

Counting runs **per sample** as its own task. Region errors are handled as described in R-12.

### 6.3 Treatment analysis (source: `Figure3/mpra_sites/modification_analysis.R` + `incell_analysis.R::diff_equiv_analysis`)
1. Subset the rows by `subset`.
2. **Master summary** per `(chr,pos)`: `avg_delrate_BS` = mean delrate of BS rows, `avg_delrate_input` = mean of input rows, `all_below_thresh = all(delrate < sesoi, na.rm=TRUE)`, `delta_delrate = BS − input`.
3. Sites with `all_below_thresh` are excluded from testing.
4. **Per tested site** (seeded `furrr::future_map`, workers = `task.cpus`; R-03):
   - `treat` is a factor with levels `c("input","BS")`.
   - If `nrow < 4`, or fewer than 2 treat levels, or fewer than 2 unique `rep`, the status is `"Not enough data"`.
   - Fit `bglmer(delrate ~ treat + (1|rep) [+ random_effects], family="binomial", weights=totalReads, fixef.prior=normal)`. On error, the status is `"Model fitting failed"`.
   - `p.value = car::Anova(model, type="III")["treat","Pr(>Chisq)"]`.
   - Set `b = mean(input delrate)`. If `b ≥ 0` and `b + sesoi < 1`, then `eqbound = logit(b+sesoi) − logit(b)`; otherwise NA. When `b = 0` the bound is infinite (D25, keep).
   - If `eqbound` is not NA: take the 90% Wald CI of `treatBS`. `is_equivalent = lower > −eqbound && upper < eqbound`. `equivalence_status` uses the detailed `case_when` labels from `modification_analysis.R`: `Equivalent by TOST`, `Different (Positive)`, `Different (Negative)`, `Inconclusive (near threshold)`, `Unmodified`, `Inconclusive (high variance / low power)`. If `eqbound` is NA, the status is `Eq. test not run`.
5. `p.adjust_diff = p.adjust(p.value, "BH")` over the tested sites only.
6. **Category** (first match wins):
   1. `all_below_thresh` → `Unmodified`
   2. `p.adjust_diff < fdr & delta_delrate > sesoi` → `Modified`
   3. `is_equivalent` → `Unmodified`
   4. otherwise → `Inconclusive`
7. **Outputs** (as in `modification_analysis.R`):
   - Full table; summary and raw tables for Modified and Unmodified sites.
   - Quartiles of Modified sites by `ntile(delta_delrate, 4)`, only if there are at least 4 Modified sites, with per-quartile tables.
   - The standardized input table and a log file.
8. **Plots** (PDF + PNG, same themes and colors as the paper code):
   - Quartile boxplot.
   - Modified-site boxplot (input vs BS averages).
   - Heatmap with and without labels.
   - All-sites faceted barplot, if enabled (the PNG may be large; that is acceptable).

### 6.4 Factor analysis (source: `incell_analysis.R::perform_anova_test` + `run_factor_dependency_analysis`)
1. Subset the rows. Make `factor` a factor with `levels` = `[baseline, experimental]`. Add indicators: `BS.ind = treat=="BS"`, `L1andBS.ind = BS & baseline`, `L2andBS.ind = BS & experimental`. There is **no low-rate pre-filter** (the legacy factor analysis has none).
2. **Per site:** if `nrow < 2`, or fewer than 2 treat levels, or fewer than 2 reps, or fewer than 2 factor levels, the note is `"Not enough data"`.
   - Full model: `bglmer(delrate ~ (1|rep) [+RE] + factor + L1andBS.ind + L2andBS.ind, ...)`.
   - Reduced model: `bglmer(delrate ~ (1|rep) [+RE] + factor + BS.ind, ...)`. Both use `family="binomial"`, `weights=totalReads` and `fixef.prior=normal`.
   - Compare them with `anova(full, reduced, test="ChiSq")` and keep the row named `glm_model_factor`. Preserve the legacy variable names so that the row name is identical.
3. `p.value.BH = p.adjust(p.value, "BH")` over the rows kept.
4. **Effect size:** take the mean delrate per `(chr,pos,level,treat)`. `delta_delrate_<level> = BS − input`, and `dd_delrate = experimental − baseline`.
5. **Significant:**
   - `p.value.BH < fdr`, and
   - either `dd_delrate > sesoi` (`positive`) or `|dd_delrate| > sesoi` (`both`).
   - For `both`, also write the split tables `<name>_<baseline>_high` (dd < 0) and `<name>_<experimental>_high` (dd > 0).
6. **Plots** (PDF + PNG):
   - Summary boxplot and heatmap of significant sites (BS minus input per level).
   - All-sites plot, if enabled.
   - Quartile boxplot when `quartiles_by` is set in a site set.

### 6.5 Site sets (source: `incell_analysis.R` union and intersection logic)
- `union` does a full join on `(chr,pos)`, with the columns suffixed `_<analysis>`.
- `intersection` does an inner join; `difference` does an anti-join.
- Optional quartiles use `ntile` of another analysis's column, as in `run_quartile_analysis`.

### 6.6 Background QC (new, D12)
Runs only when `count_all_bases` is set.

- **Per sample and ref base:** the distribution of delrate, reported as median, IQR and mean.
- **Per ref base:** input vs BS.
- **Outputs:** a table, plus a boxplot faceted by ref base (PDF + PNG).

Site calling always filters `ref == "T"`, so this mode never changes the calls.

---

## 7. Work packages

Dependency graph: **WP0 → (WP1 ∥ WP2 ∥ WP3 ∥ WP4 ∥ WP5 ∥ WP6) → WP7 → WP8.**

The Orchestrator runs WP0 itself. It then spawns one agent per WP in WP1–WP6, each in an isolated git worktree. It merges PRs, keeps `PROGRESS.md` current, and runs WP7–WP8.

### WP0 — Scaffold and contracts (Orchestrator)
- Create the nf-core pipeline template (`nf-core pipelines create`, name `nanomodamp`) inside the NanoModAmp repo. Keep the LICENSE.
- Write `nextflow_schema.json` (§3), `assets/schema_input.json` (§5.1), `assets/schema_analyses.json` (§5.4), and `docs/contracts/*.md` (§5.2–5.5, with example files).
- Stub every process so that `nextflow run . -profile test,docker -stub` completes end to end.
- Set up CI (GitHub Actions): nf-core lint, nf-test, R `testthat`, pytest, and the `-profile test` run on Docker. The Singularity run goes in a separate optional job.
- Create `PROGRESS.md`, `CHANGE_REGISTER.md` (seeded from §9), `docs/DECISIONS.md` (seeded from §1) and `CHANGELOG.md`.
- **Accept when:** the stub run passes and CI is green. **[GATE G0]** Becca reviews the contracts and schemas.

### WP1 — Synthetic simulator and truth data (Python)
- `tests/simulate/simulate.py` takes a small reference (5 amplicons + 20 MPRA-like oligos, all `+` strand) and a BED, and generates per-sample FASTQs with:
  - Both orientations, ONT flanks, 10-nt 3′ UMIs, and MPRA pool adapters when `--library mpra`.
  - Specified per-site deletion probabilities for each sample (from a YAML design).
  - Replicates, and effects for BS vs input and for WT vs KD.
  - PCR duplicates at a known rate, including 1-mismatch UMI copies.
  - Junk reads: no adapter; **only one adapter**; too short; adapter dimers; chimeras with both orientations; multi-mappers across two near-identical oligos; insertions next to a T; low-coverage sites (≤ 20 reads); sites where input delrate = 0.
  - An optional error model (Badread), switched off by default so tests are deterministic.
- Output: FASTQs, sample sheet, analyses YAML, and **truth tables** (expected funnel counts per filter, expected per-site counts after dedup, and expected categories for clear-cut sites).
- Commit the generated data to `tests/data/synthetic/` (≤ 5 MB), along with the generation command and seed.
- **Accept when:** the generation is deterministic (same seed gives identical bytes), pytest checks pass on the simulator itself, and the design is documented.

### WP2 — Preprocessing subworkflow (Nextflow)
- Use nf-core modules where available: cutadapt, seqtk/seq, umitools/extract, minimap2/align, samtools, umicollapse. Put arguments in `conf/modules.config` and write local modules for the rest (CAT_FASTQ, MERGE_ORIENT, READ_FUNNEL).
- Pin all container versions and emit `versions.yml`.
- **nf-test cases:**
  1. A sense read is trimmed to exactly the insert.
  2. An antisense read comes out identical to its sense twin after RC.
  3. **A read with only one adapter is discarded** (D4).
  4. A too-short read is discarded.
  5. The extracted UMI equals the simulated UMI, in the read-name format UMICollapse expects.
  6. Duplicates and 1-edit UMIs collapse; distinct UMIs are kept.
  7. Secondary, supplementary and MAPQ < 30 reads are removed, and MAPQ = 30 is kept.
  8. A directory input concatenates deterministically.
  9. An empty directory or missing file gives a clear error.
  10. The raw `.gz` count is correct.
  11. Every read-funnel number matches the WP1 truth.
  12. `-profile test` works with both library types.
  13. `orientation_adapters=pool, umi=false` (D29): a pool-flanked read without ONT adapters or UMI is kept, oriented and mapped; no UMI step or dedup runs; the funnel omits the skipped steps; `pool` + `umi=true` and `pool` with `library_type=endogenous` fail validation.
  14. `ont_adapter_mode=three_prime_only` (D31): a read with only the 3′ ONT adapter is kept and trimmed at the 3′ end; output equals legacy `cutadapt -m 125 -O 15 -a GAAGATAGAGCGACAGGCAAGT` on the same reads; the default `linked` discards it.
  15. (WP3) `input_counts` (D32): two tables merge into the same result as one table with all samples; duplicate `sample_id` and mismatched columns fail clearly; shuffled column order gives the same result.
- **Accept when:** all tests pass, and the minimap2 defaults resolve by library type and can be overridden.

### WP3 — Counting, merge and background QC (R)
- Implement §6.2 in `nanomodamp::count_sites()`, with the CLI `bin/nma_count.R`. Then write `merge_counts()` (§5.3, joining by name) and `background_qc()` (§6.6).
- Region errors are handled as described in R-12, `max_depth` saturation gives a warning (R-16), and the `bed_coordinates` parameter is supported (R-15).
- **testthat cases (using WP1 BAMs, or tiny BAMs built with Rsamtools in a fixture):**
  1. The deletion count equals the truth.
  2. T-only filtering works.
  3. All-bases mode works.
  4. kmer is NA at the edges.
  5. A `-` strand region is complemented and reverse-complemented correctly (defensive, even though D19 says strand is always `+`).
  6. Wrong-strand reads are excluded.
  7. For single-site BEDs vs amplicon BEDs, positions are correct under both `bed_coordinates` options.
  8. An insertion next to a T adds exactly 1 to that position's `totalReads` (this documents D9).
  9. A site with exactly 20 reads is excluded and one with 21 is kept.
  10. Merging tables with shuffled column order gives the same result.
  11. Rerunning produces identical output (idempotent).
  12. A thread count of 1 vs N gives identical output.
- **Accept when:** all tests pass, and on WP1 BAMs the per-site counts equal the truth exactly.

### WP4 — Site calling, site sets and plots (R)
- Implement §6.3–6.5 as `treatment_test()`, `factor_test()`, `classify_sites()`, `site_sets()` and the plot functions, with the CLIs `bin/nma_call.R` and `bin/nma_sitesets.R`. Every plot is saved as PDF + PNG.
- Build `assets/analyses_example.yaml`, which reproduces the paper analyses (Figure 3 in vitro, the Figure 3 in-cellulo WT_mod / PUS7_dep / celltype_spec runs for HepG2, 293T and Both, and the Figure 2 endogenous runs), translated into generic treatment/factor entries.
- **testthat cases:**
  1. A clear positive site is `Modified`, a clear null is `Unmodified`, and a low-coverage noisy site is `Inconclusive`.
  2. The all-below-SESOI pre-filter works.
  3. `Not enough data` and `Model fitting failed` paths are handled.
  4. `p.adjust` receives exactly one p-value per tested site.
  5. A site with input = 0 behaves as the legacy code does (D25).
  6. A simulated interaction is detected for `positive` and `both`, and the sign split is correct.
  7. A site with equal effects is not significant.
  8. Seeded runs on 1 core and on N cores are identical.
  9. The schemas match §5.5.
  10. Every expected plot file exists, is non-empty, and has both PDF and PNG versions.
  11. **Exact equivalence with legacy:** given the same input table, `treatment_test()` matches `modification_analysis.R`, and `factor_test()` matches `incell_analysis.R`, in p-values (tolerance 1e-8), categories and dd values. Use the paper's shipped `Figure3/parameter_sweep/BIDdetect_data_invitro_delpos.txt` and `Figure3/plots/BIDdetect_data_incell_delpos.txt`, and compare against the shipped `invitro_modification_significance.tsv`, `WT_mod_Both_significance.tsv` and `PUS7_dep_union_significant_summary.tsv`. This test is independent of WP5.
- **Accept when:** all tests pass, including test 11 with zero unexplained differences.

### WP5 — Legacy harness and golden test package (Python/R; needs Becca)
**Detailed plan:** `golden_test_package_plan.md`. Decisions are logged in `docs/DECISIONS.md`. Updated 2026-10-06 to reflect the G1 decisions below.

1. **Two harnesses.**
   - `dev/legacy/`: the scripts **actually executed** for each paper run (G1-a/F1). The paper numbers were not produced by the published Figure 2/3 scripts. They came from earlier versions using minimap2 `-a -k5`, no primary/MAPQ filter, and `umi_tools dedup --method directional`. In cellulo had no UMI and no dedup. Allowed edits: (a) keep the pre-dedup BAM, (b) paths/modules/no runtime installs, (c) in-cellulo config, (d) `sample_name.R` (D26), (e) a fixed `umi_tools dedup --random-seed` (G1-c; the paper runs set none, so their dedup representatives cannot be regenerated). `ORIGINAL_DIFF.patch`, `SOURCES.md` and `versions.txt` record everything.
   - `dev/published/`: the published PUS7regulation2026 Figure 2/3 scripts (UMICollapse, `-ax sr` / `-ax splice -uf`, `-F 2304 -q 30`), same allowed edits. In cellulo uses a documented constructed variant (pool-adapter orientation + published mapping/filters, no UMI), because the published scripts have no in-cellulo path (D29, R-26).
   - Containers for both harnesses (`containers/legacy/`) follow once the Sherlock runs are frozen; tool versions are in `versions.txt`.
2. **Targets (G1-b):** endogenous RHBDD2, HDAC6 (run 20250418) and STIM1 (run 20251022) with `endoG1G2.fasta` (F2); MPRA: the six oligos in `pool1_chr.txt` for in vitro and in cellulo. Accepted gaps: no Inconclusive, input-delrate = 0, or Modified-not-PUS7-dependent MPRA site. Cell-type-specific analyses are out of scope; tests focus on treatment and PUS7 dependency.
3. **Samples (G1-b):** endogenous 32 (20250418 reps 1, 3; 20251022 reps 1, 2; P36/pLKO → WT, P97/shPUS7 → KD), in vitro 4 (May run; noPUS → `treat=input`, D14), in cellulo 32 (HepG2 reps 1, 2; 293T reps 1, 3; P101/P102/P3/P4 × input/BS).
4. **Reads:** all reads with any alignment on a selected target in the rerun pre-dedup BAM, with their duplicates, extracted raw from the original FASTQs in original order. **Size (G1-d):** at most **250 selected reads per (sample, target)** (seeded, nested) plus **150 background reads per sample** (seeded); total ≤ 25 MB. Targets below the cap are kept whole.
5. **Checks:** L1 (rerun BAM vs paper BAM), L2 (extraction complete), V1 (legacy counts on the full selection vs paper counts: **exact** for runs counted on non-deduplicated BAMs — 20250418, in cellulo; **within tolerance** for runs counted on dedup BAMs — 20251022, in vitro: same `totalReads`, `*.count` within max(3, 5%), delrate within 0.03), Repro (byte-identical rebuild).
6. **Expected outputs** under `tests/data/golden/<experiment>/expected/`:
   - `paper_reference/` (reference only), `legacy_full_selected/` (V1 counts),
   - `legacy_rerun/` (**binding** for counts; endogenous calls reference only, glmer),
   - `published_rerun/` (D30; reference for the WP7 comparison).
7. **R-15 (BED convention):** the delpos BEDs used for the paper (`endoG1G2_delpos.bed`, pool1 delpos) have start = end = the 1-based site, and the legacy counter adds no +1, so counts are correct for them. Standard 0-based BEDs (e.g. `set3_delpos.bed`) would be off by one under the legacy code. Evidence and recommendation in `dev/golden/R15_bed_convention.md` **[GATE G1-e]**.
- **Accept when:** L1, L2, V1 (per the rule above), V2 and Repro pass; both reruns are stored; Becca approves G1-a through G1-e; the package is committed via PR to `dev`.

### WP6 — Documentation
Write `docs/usage.md` (quick start with `-profile test,docker`, a sample sheet guide, the analyses YAML guide, and running on SLURM, Sherlock and cloud), `docs/parameters.md` (generated from the schema), `docs/output.md` (every file and column), and `docs/methods.md`. `methods.md` must cover:

- **Why `-ax sr` is the default for MPRA and `splice -uf` for endogenous** (D5): these are validated in the paper, alignment presets change how deletions are represented, and the page must explain how and when to change them.
- The denominator note (D9).
- Coverage filtering (D10), including the legacy parameter-sweep's `>= 20` wording.
- The BED convention (R-15).
- The treatment/TOST/factor models, written out as equations.
- The equivalence limitation at input = 0 (D25).
- The replicate pairing assumption (D17).
- Requiring both adapters (D4).
- Libraries without ONT adapters or UMIs (D29): state that the defaults (ONT + UMI + dedup) apply to all library types including in cellulo; the opt-in `orientation_adapters=pool, umi=false` exists only for legacy libraries like the paper's 20241114 in-cellulo run (used by `test_golden`); what changes in the funnel and outputs. Same for `ont_adapter_mode=three_prime_only` (D31): opt-in only for libraries lacking the 5′ ONT adapter, like paper run 20250418.
- How the golden test works.

Also write `CITATION.cff`, and a README with a badge for the test profile.

- **Accept when:** a new user can follow `usage.md` from a clean machine with Docker and reproduce the `test` and `test_golden` results.

### WP7 — Integration and golden concordance (Orchestrator)
1. Merge WP1–WP6 and run `-profile test` and `-profile test_golden` on Docker and Singularity.
2. **Concordance report** (`docs/validation/golden_concordance.md`, generated by `dev/concordance.R`): compare stage by stage against `expected/legacy_rerun`, `expected/published_rerun` (D30) and `expected/paper_reference`. Differences vs `published_rerun` isolate the new implementation; differences between `legacy_rerun` and `published_rerun` isolate the paper-processing changes (R-25). Cover:
   - Read funnel.
   - Per-site `totalReads` and `delrate`: exact-match rate, max |Δ|, and a scatter plot.
   - Category confusion matrices for each analysis.
   - Factor significance agreement.
3. **Attribute every difference** to a Change Register entry (for example R-04 requiring both adapters, or the modern cutadapt version, or D13 bglmer for endogenous, or R-15). Rerun with `require_both_adapters=false` and the legacy tool versions to confirm each attribution. Any difference that can't be attributed blocks the release.
4. **[GATE G2]** Becca reviews the report. Once approved, freeze the pipeline outputs as `tests/data/golden/expected/pipeline_v1/`. Add an nf-test that checks exact regression against the frozen outputs, with numeric tolerance 1e-8 and exact categories.
- **Accept when:** G2 is approved and CI runs the golden regression test.

### WP8 — Release
Version `1.0.0`. Finalize `CHANGELOG.md` and the Change Register. Make a Zenodo release (optional, for a DOI). Pin the container digests. Do a final pass with `nf-core pipelines lint`. **[GATE G3]**

---

## 8. Test matrix summary

| Layer | Tool | Data | Where |
|---|---|---|---|
| Simulator | pytest | — | WP1 |
| Modules/subworkflows | nf-test | synthetic | WP2, WP3 |
| R functions | testthat | fixtures + synthetic | WP3, WP4 |
| Stats equivalence vs paper code | testthat | shipped paper tables | WP4 |
| End-to-end | nf-test pipeline, `-profile test` | synthetic, both library types | CI |
| Golden QC for users | `-profile test_golden` | paper subset | CI + users |
| Golden regression | nf-test snapshot | frozen `pipeline_v1` | CI after G2 |
| Portability | CI matrix | Docker (required), Singularity (optional job), conda (smoke) | CI |

---

## 9. Change Register (seed for `CHANGE_REGISTER.md`)

Status values: `APPROVED` (implement), `KEEP` (preserve legacy behavior and document it), `SUPERSEDED` (not applicable in the new design), `OPEN` (decide at a gate).

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
| R-28 | `umi_tools dedup` picks random representatives (no seed) | Superseded by UMICollapse (D7); legacy harness seeded (edit e) | SUPERSEDED | None |

---

## 10. Gates and open questions

| Gate | When | Becca decides / provides |
|---|---|---|
| G0 | After WP0 | Approves the contracts (§5), parameter names and defaults (§3), and the analyses YAML format |
| G1 | During WP5 | Provides local data paths and exact tool versions used in the paper (minimap2, UMICollapse commit, umi_tools, seqtk, R/package versions); approves target and sample selection; decides R-15 |
| G2 | After WP7 | Approves the concordance report; authorizes freezing `pipeline_v1` golden outputs |
| G3 | Release | Approves v1.0.0 |

Open items to confirm at G1:

- The exact legacy tool versions.
- The BED convention (R-15).
- Whether the endogenous analyses in `analyses_example.yaml` should use `(1|vector)` or `(1|celltype)` random effects, to mirror the Figure 2 design as closely as possible under the Figure 3 method.
- Whether reads found in both orientations exceed 0.1% of reads.

Resolved at G1 on 2026-10-06 (see `docs/DECISIONS.md`): tool versions (`dev/legacy/versions.txt`), targets, samples, size policy, harness scripts (F1), endogenous reference (F2), dedup seed and V1 rule (G1-c). R-15 decided at G1-e the same day. Endogenous random effects decided 2026-10-09 (D33).
