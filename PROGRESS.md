# Progress

## WP5 — Golden test package (`wp/5-golden-package`)

| Step | Status |
|---|---|
| G1-a preflight findings (`NanoModAmp_build/preflight_findings.md`) | Done; F1/F2/F4 resolved 2026-10-06 |
| Preflight checks 0–6, input provenance | Done — reports in `NanoModAmp_build/golden/reports/` |
| §4.1 legacy harness (`dev/legacy/`) | Done (not yet run); see `dev/legacy/SOURCES.md` |
| G1-b chr coverage profile + minimal sample list | Approved 2026-10-06 |
| G1-d size/depth | Approved 2026-10-06: cap 250 selected reads per (sample, chr) + 150 background per sample |
| §4.2 L1 legacy rerun on approved samples | Submitted 2026-10-06: arrays 46826528/30/33/38, compare 46826569 (`dev/golden/submit_l1.sh`) |
| §4.3–4.5 selection, background, extraction (cap 250 / bg 150) | Submitted: job 46827301 after L1 (`dev/golden/run_select.sh`, `select_extract.py`) |
| §4.6 V1 legacy counts on full selection vs paper | Submitted: job 46828547 after selection (`dev/golden/run_v1.sh`, `v1_compare.py`) |

### Open questions
- Endogenous vector labels normalised in sample sheets: P36/pLKO → WT, P97/shPUS7 → KD (Becca 2026-10-06); original label kept in `notes`.
- 20250418 paper counts are pre-dedup (`sorted/`, MAPQ≥1 in counting); 20251022 are dedup. README must state the mixed source.

### L1 test findings (sample 293T_P36_1_BS, 2026-10-06)
- Harness reproduces the paper pre-dedup BAM **byte-for-byte** (alignments + sequences) with the original reference, at 4 or 15 threads.
- With endoG1G2.fasta (F2): same read IDs; some supplementary/CIGAR records differ; pileups at RHBDD2:287 and HDAC6:392 unchanged.
- umi_tools dedup is non-deterministic (paper runs set no `--random-seed`): same number of reads kept, different representatives; site base counts shift by a few reads. Affects counts of 20251022 and in vitro (counted on dedup BAMs); not 20250418 (counted pre-dedup) or in cellulo (no dedup).
- G1-c resolved 2026-10-06: harness edit (e) `--random-seed` added (verified deterministic); V1 tolerance rule for dedup-counted runs (see DECISIONS.md).
- Selection/extraction tested: L2 pass, byte-identical on rerun, ~294 B/read gz; rerunning on the extracted reads gives alignments identical to the full run.
- **Stale index (2026-10-06):** `/home/groups/nicolemm/rodell/pool1/pool1_cleaned_noadapters.fasta.fai` (2025-11-07) predates the FASTA's CRLF→LF rewrite (2026-01-26); offsets are wrong (e.g. line 1 offset 26 vs 25). Golden drivers stage a copy with a fresh `.fai` (`dev/golden/stage_refs.sh`); the original is untouched. Becca to decide whether to rebuild it in place (`samtools faidx`). `endoG1G2.fasta.fai` verified current.

## Resume notes (written 2026-10-06 16:35; interactive allocation 46816746 ends 22:13)
Slurm jobs run independently of the interactive session. Chain and outputs (all under `/scratch/users/rodell/NanoModAmp_build/golden/`):

| Job | What | Depends on | Output / check |
|---|---|---|---|
| 46826528, 46826530, 46826533, 46826538 | L1 legacy reruns (arrays) | — | `l1/<exp>/`, logs in `l1/slurm_logs/` |
| 46826569 | L1 compare | afterok L1 | `l1/L1_results.tsv` |
| 46827301 | Selection / extraction | afterok L1 | `select/select_summary.tsv`, `select/select_*.out` (sizes) |
| 46828547 | V1 | afterok select | `v1/V1_results.tsv`, `v1/v1_*.out` (PASS/FAIL) |
| 46830631 | V2 legacy rerun + size check | afterok V1 | `expected/legacy_rerun/`, `size_report.txt` |
| 46831003 | Published rerun | afterok select | `expected/published_rerun/<exp>/status.tsv` |
| 46837506 | Repro (selection rebuild, byte-identical) | afterok select | `select/repro_*.out` (PASS/FAIL) |
| 46839079 | Package assembly → `tests/data/golden/` (+ size check, MANIFEST.md5) | afterany V2, published | `dev/golden/run_assemble.sh` log; partial build in `assemble_partial/` if inputs missing |

If a job fails, dependents stay PENDING with reason `DependencyNeverSatisfied`; find the cause (`sacct -j <id>`, `.err` logs) before resubmitting.
Next steps after the jobs finish: review L1/V1 (calibrate the V1 tolerance), check size ≤ 25 MB, §4.7 package assembly (sample sheets with WT/KD and run-unique reps as in `expected/published_rerun/endogenous/name_map.tsv`), §4.10 manifest + Repro rebuild, PR to `dev`.
WP0: worktree `/home/users/rodell/NanoModAmp_wp0`, branch `wp/0-scaffold` (pushed to GitHub).

- G1-e decided 2026-10-06 (bed0 default; golden targets.bed 0-based). Dry-run leftovers safe to delete: `golden/assemble_dryrun/`, `$GROUP_HOME/rodell/envs/nfcore` (broken py3.6 venv).

### Results 2026-10-06 18:25
- **L1** (`golden/l1/L1_results.tsv`): pre-dedup read IDs identical to paper for every sample × target (32/32, 16/16, 24/24; in cellulo 192/192 incl. alignments). In vitro dedup BAMs identical too (24/24). 20251022 pre-dedup alignments identical; dedup keeps the same read count (16/16) with different representatives (expected, unseeded paper umi_tools). 20250418: alignments differ for some reads because of endoG1G2.fasta (F2); dedup counts differ by ≤ 4 reads in 10/32 (positions shift) — no effect on its counts, which are pre-dedup.
- **L2:** pass for all 68 samples. **Repro (selection):** PASS, 340 files byte-identical.
- **Size:** capped FASTQs 19.2 MB (endo 0418 3.31, 1022 1.90, in cellulo 12.03, in vitro 1.98); full selection 72.8 MB (V1 input only).

### Resubmission 2026-10-06 18:55
First V1 (46828547) and published (46831003) runs failed: Slurm snapshots the batch script at submission, so both ran versions from before later fixes (V1 lacked `stage_refs.sh` → stale pool1 .fai → in vitro counting failed; published lacked the D31 variant and had an R 4.3.2 module purge that dropped libgfortran). V2 (46830631) was cancelled by dependency; assemble (46839079) failed on missing expected outputs (size without them 19.57 MB). Failed outputs moved to `golden/failed_*`. Fixed, tested (published in vitro site calling OK), and resubmitted:

| Job | What | Depends on |
|---|---|---|
| 46844253 | V1 (`v1/V1_results.tsv`) | — |
| 46844259 | V2 legacy rerun | afterok V1 |
| 46844262 | Published rerun | — |
| 46844265 | Package assembly | afterany V2, published |

### V1 accepted 2026-10-06 19:50 (G1-c); resubmitted V2 46852559 and assembly 46852562 (afterok V2). Published rerun 46844262 completed ok.

### Golden package built 2026-10-06 20:35
`tests/data/golden/` assembled (job 46852562): 251 files, **22.3 MB** (cap 25 MB; `du` shows ~49 MB only because of Isilon block allocation).
Checks: L1 ✓, L2 ✓ (68/68), V1 ✓ (242/242 rows, refined rule), V2 ✓ (legacy rerun 46852559, size check PASS), published rerun ✓ (46844262), Repro ✓ (selection 340 files; package MANIFEST.md5 byte-identical on rebuild, 250 files). Gates G1-a…G1-e approved.
Legacy vs published site calls identical for MPRA: in vitro STIM1/RHBDD2/HDAC6 Modified, ZNF644/TRIM71 Unmodified, GOSR1 Inconclusive; in cellulo WT_mod_Both STIM1/RHBDD2/HDAC6 Modified; PUS7_dep union GOSR1/HDAC6/RHBDD2/STIM1.
Remaining for WP5 done: commit the branch and open the PR to `dev` (needs a `dev` branch; WP0 merge first). Endogenous random effects decided 2026-10-09 (D33).
