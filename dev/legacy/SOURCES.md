# Legacy harness: script sources

Under Becca's F1 decision (2026-10-06), the harness reruns the scripts that were **actually executed** for the paper runs, not the published PUS7regulation2026 Figure 2/3 `file_prep/` scripts. `original/` holds pristine copies. Working copies in `endogenous/`, `mpra_invitro/` and `mpra_incell/` differ from them only by the edits in `ORIGINAL_DIFF.patch`. Source `env.sh` before running anything.

Repos used:
- `rmrodell/BIDamplicon`, local clone at `/home/users/rodell/BIDamplicon`, HEAD `315d840` (2026-04-04). The clone has uncommitted edits in `pipeline/` (trim_map_dedup*.sh, incell_analysis.R) and `legacy/DelDetect.sh`. None of these was used: every file below came from `git show <commit>:<path>` or from a run directory.
- `martinezlab/PUS7regulation2026`, local clone at `/home/users/rodell/PUS7regulation2026`, HEAD `4fd285b` (2026-04-24). Site-calling scripts were last changed in `4673e50` (2026-04-17), and the clone is clean.

**How a version was matched.** Each run called the scripts live from `$HOME/BIDamplicon/pipeline/` (the working tree). For each script, the version taken is the latest commit at or before the job start time. That choice was then checked against command lines and log formats in the run logs. Where the working tree was ahead of git at run time, the next commit is used and the reason is stated.

## Endogenous run 20250418 (processing `legacy_dedup_mapping`)

The paper data come from the rerun in `.../20250418_endoBID/20250503_70/`, not from `prep_map.sbatch` / job 63636074. That earlier job wrote `.../20250418_endoBID/deduplication/`; the yaml points at `20250503_70`.

| File (original/endogenous/20250418/) | Source | md5 | mtime | Executed by | Evidence |
|---|---|---|---|---|---|
| 20250503.sbatch | run dir `20250503_70/` | aa762cc1… | 2025-05-10 | job 64514787 (2025-05-03): fastaprep, mapping_4psi, dedup_map_folder, DelDetect | `.out`/`.err` in run dir |
| 20250503_full.sbatch | run dir | a343be69… | 2025-05-10 | job 64526923: DelDetect on `deduplication/minimap2/sorted` with `20250404_endoBIDamplicon_full.bed` | `.out` ends "Final count file located at …/20250503_70/deduplication/counts/deldetect_counts.txt" |
| fastaprep.sh | BIDamplicon `d77c723:pipeline/fastaprep.sh` (unchanged through 181d080) | 83e5df19… | — | job 64514787 | log "MD5 checksum for X.fq matches". **Not run by the harness**: it runs `gunzip` inside the source directory. The harness starts from its output, `20250503_70/fastq/*.fq`, which is the yaml `fastq_dir`. |
| dedup_map_folder.sh | `d77c723:pipeline/dedup/dedup_map_folder.sh` | 0e4f5c73… | — | job 64514787 | identical at d77c723, 181d080 and HEAD legacy/ |
| dedup_mapping.sh | `d77c723:pipeline/dedup/dedup_mapping.sh` | 5edf81b8… | — | job 64514787 (2025-05-03) | logs: `cutadapt 1.18 -m 125 -O 15 -a GAAGATAGAGCGACAGGCAAGT`; `umi_tools extract --extract-method=string --bc-pattern=NNNNNNNNNN --3prime`; `minimap2 -a -k5 -t 15`; `umi_tools dedup --method directional`. The next commit (181d080, 2025-05-06) only adds `-F 4` to the summary read counts, so BAMs are identical under either version. |
| DelDetect.sh | `d77c723:pipeline/DelDetect.sh` (unchanged until 2025-07-23) | 34c47387… | — | job 64526923 (~2025-05-10) | `.out` "Processing sample: 293T_P36_1_BS_UMI" |
| bam_counts.R | `d77c723:pipeline/bam_counts.R` (unchanged until 2025-07-23) | 7328ec82… | — | job 64526923 | serial version; `.out` "Examining HDAC6 , number 3 out of 4 total genes." |
| sample_name.R | `d77c723:pipeline/sample_name.R` (= HEAD `legacy/sample_name.R`, D26) | 48f9939d… | — | job 64526923 | `deldetect_factors.txt` has the celltype/vector/rep/treat split |

**Important for L1/V1:** the 20250418 paper counts came from the **pre-dedup** `sorted/*_UMI_sort.bam`, not from `dedup/`.
- The sbatch passes `.../minimap2/sorted` as the BAM directory.
- The sample column is `293T_P36_1_BS_UMI`. That is `basename X_UMI_sort.bam _sort.bam`; a dedup BAM would have produced `..._UMI_dedup.bam`.
- The paper table `20250418_BIDdetect_data.txt` matches the run's `deldetect_factors.txt` row for row: all 156 non-4psi rows are identical (HDAC6:390 293T_P36_1_BS gives totalReads 78 in both). Its 92 `4psi` rows came from the separate `4psi.fasta` mapping, which is out of scope.
- Counting used `PileupParam(min_mapq = 1, min_base_quality = 1)`, so MAPQ-0 reads are dropped at counting. That likely explains 151 sorted reads becoming 78 counted.

## Endogenous run 20251022 (processing `fig2_trim_map_dedup`, run-dir copy)

| File (original/endogenous/20251022/) | Source | md5 | mtime | Executed by | Evidence |
|---|---|---|---|---|---|
| prep_fastq.sh | run dir `20251022_endoBID/` | 1f3e73a4… | 2025-12-18 (copy to Oak) | interactive, per process.txt Step 1 | output `fastq/<sample>.fastq.gz` is on Oak; concatenation is `cat barcodeNN/*.fastq.gz` (glob, lexical order) |
| submit_trim_map_dedup.sh | run dir | a2d77a37… | 2025-12-18 | submitted job 8623868 | process.txt transcript "Submitted batch job 8623868", job name SHAPE_prep |
| trim_map_dedup.sh | run dir; **identical to BIDamplicon `4b3900f`** (committed 2025-10-30, the day after the run) | fea5429f… | 2025-12-18 | job 8623868 (2025-10-29 12:44) | `mapping/*/logs/*.pipeline_run.log` "[Job:8623868 Task:10] …"; BAM @PG `minimap2 -a -k5 -t 16 …/set3.fa` |
| BIDdetect.sh | `9c6e9ba:pipeline/BIDdetect.sh` (committed 2025-10-29 14:44) | c6aa0392… | — | job 8640004 (started 2025-10-29 14:49:33) | `.out` prints "Final Column Names:" (added in 9c6e9ba) and calls `bam_counts.R` (the switch to `_fast` came 2025-11-06) |
| bam_counts.R | `9c6e9ba` (= f500189, parallel + delrate) | bbd08f19… | — | job 8640004 | `.out` "Starting counting at 2025-10-29 14:49:44" (only in the f500189+ version) |
| sample_name.R | `9c6e9ba:pipeline/sample_name.R` | 48f9939d… | — | job 8640004 | — |
| count.sbatch, process.txt | run dir | b0b5e3d0…, e324ddb2… | 2025-12-18 | record only | count.sbatch on Oak is the later `counts_full` edit; the delpos call is recorded in process.txt |

Paper table `endoBID_data/20251022_BIDdetect_data.txt` is md5-identical (06b72d48…) to `20251022_endoBID/counts_delpos/BIDdetect_data.txt`. The pre-dedup BAM `mapping/<s>/tmp/<s>_mapped_sorted.bam` is still on Oak.

## MPRA in vitro, runs may (job 9257951) and june (job 9257956)

| File (original/mpra_invitro/) | Source | md5 | Executed by | Evidence |
|---|---|---|---|---|
| submit_file_prep.sh | `ab9db01:pipeline/submit_file_prep.sh` (only version) | 7cc0639a… | submit of 9257951 / 9257956 | process.txt transcript |
| trim_map_dedup_mpra.sh | `c3ca8e4:pipeline/trim_map_dedup_mpra.sh` | fa91e19d… | 9257951 (2025-11-07 11:49), 9257956 | The job ran before c3ca8e4 was committed (13:34), but its `slurm_logs/*.out` print "Info: Found barcode subdirectory … Concatenating multiple FASTQ files into …/00_concatenated_fastqs/…", which exists only from c3ca8e4 on; ab9db01 would have looked for `<sample>.fastq.gz` and failed. So the working tree already held the c3ca8e4 content. Cutadapt command lines in logs (`-m 125 -O 15 --discard-untrimmed -g TTTCTG…...GAAGAT…`; pool trim `-m 120 -O 10 -g GACGCTCTTCCGATCT...CACTCGGGCACCAAGGAC`) match. |
| BIDdetect.sh | `2046eb2:pipeline/BIDdetect.sh` (2025-11-06 19:44) | 2ff62568… | job 9275304 (2025-11-07 14:26), which wrote `counts_delpos/logs/BIDdetect_20251107_142611.log` | calls `bam_counts_fast.R`; the next change is 2025-11-26 |
| bam_counts_fast.R | `bf170b3` (2025-11-06 14:43) | 9f29c726… | job 9275304 | `39f2fc1` (kmer edge fix) was committed 16:51, after the count job. It only changes `kmer` for pos < 3 or pos > len−2; none occur in the delpos table (0 empty kmers), so the output is the same under either version. |
| sample_name.R | `2046eb2` | 48f9939d… | 9275304 | `--col_names celltype_vector_treat_rep` |
| modification_analysis.R | PUS7regulation2026 `Figure3/mpra_sites/` (`4673e50`) = BIDamplicon `678c93a` (2025-11-17) | 8abf253f… | process.txt Step 3 | byte-identical to BIDamplicon 678c93a |
| process.txt, count.sbatch | `Pool1/invitro/` | 48c605f6…, 66f8ac7e… | record | count.sbatch on Oak is the later counts_full edit |

The concatenated input FASTQs are kept on Oak (`maydata|junedata/00_concatenated_fastqs/`). Concatenation used `zcat $(find barcodeNN -name '*.fastq.gz')`, i.e. **`find` (directory) order, not lexical**. Golden §4.5 must reuse these files, or reproduce this order, rather than re-sorting. The pre-dedup `tmp/<s>_mapped_sorted.bam` is on Oak.

## MPRA in cellulo, run 20241114 (job 9256680)

| File (original/mpra_incell/) | Source | md5 | Executed by | Evidence |
|---|---|---|---|---|
| submit_file_prep.sh | `ab9db01` | 7cc0639a… | submit of 9256680 | process.txt |
| trim_map_mpra.sh | `12581ec:pipeline/trim_map_mpra.sh` ("custom built for nov2024 pool1", 2025-11-06 20:02; unchanged in c3ca8e4) | 7e48ea76… | 9256680 (2025-11-07 11:27) | `.out` lacks the "Sample Map File:" line (removed in 12581ec) and reads `<sample>.fq` by sorted index; cutadapt logs show pool adapters as the linked pair (`-g GACGCTCTTCCGATCT...CACTCGGGCACCAAGGAC`), no UMI step, no dedup |
| BIDdetect.sh / bam_counts_fast.R / sample_name.R | `2046eb2` / `bf170b3` / `2046eb2` | as in vitro | job 9276805 (2025-11-07 14:36:27), log `BIDdetect_20251107_143627.log` | default `--col_names celltype_vector_rep_treat` |
| incell_analysis.R | PUS7regulation2026 `Figure3/mpra_sites/` (`4673e50`) | dd767531… | — | differs from BIDamplicon `51dc8f4` by 13 diff lines (config). The paper repo copy is the published one. |
| process.txt, count.sbatch | `Pool1/incell/` | 7904101f…, 99c5c6a1… | record | — |

`prep/final_bam/` has mtime 2026-03-30 because `read_counts.tsv` was added then. The BAMs are dated 2025-12-18 (Oak copy).

## Endogenous site calling

| File | Source | md5 |
|---|---|---|
| original/endogenous/site_calling/analysis_endo.R | PUS7regulation2026 `Figure2/endo_sites/` (`4673e50`) | 6145317f… |

`original/common/sample_name.R` is BIDamplicon HEAD `legacy/sample_name.R` (D26). It is byte-identical to every per-run copy.
