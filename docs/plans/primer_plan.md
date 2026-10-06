# NanoModAmp Amplicon & Primer Design Tool — Plan for Claude Agents

**Goal:** Turn the manual endogenous-amplicon design process in `rmrodell/BIDamplicon/primers/` into a reproducible standalone tool, **`nma-design`**. It starts from a list of genomic target sites. It outputs ready-to-order primer oligos plus the amplicon FASTA and BED files that the NanoModAmp endogenous workflow uses as its reference.

- **Source (read-only):** `https://github.com/rmrodell/BIDamplicon/tree/main/primers` (`amplicon_design`, `How_To_amplicon_design.txt`)
- **Target:** the `design/` directory of `rmrodell/NanoModAmp`. This is a separate work package, **WP9**, and runs in parallel with WP1–WP6 of `plan.md`.
- **Owner / reviewer:** Becca (rodell@stanford.edu). Steps marked **[GATE]** need her approval.

The agent rules from `plan.md` §0 apply here: frozen behavior unless registered, tests before merge, PROGRESS/DECISIONS logging, and stopping when unsure.

---

## 1. What the current process does

1. **`process_coordinates(coords)`** takes GRCh38 `chrN:pos` sites. For each site:
   1. Pick the strand whose base is T.
   2. Query live Ensembl (biomaRt, no pinned release) for canonical transcripts on that strand within ±1 kb that have a RefSeq mRNA ID.
   3. Take the first one.
   4. Compute the 1-based transcript position from the Ensembl exon coordinates.
   5. Set `Fprimer_to = tp − 50` and `Rprimer_from = tp + 50`.
   6. Write an Excel file.
2. **Manual PrimerBLAST**, done per site, on the RefSeq mRNA:
   - Forward primer must end at or before `Fprimer_to`, and the reverse primer must start at or after `Rprimer_from`.
   - Product 200–1000 bp; primer Tm max 72 °C, optimum in the upper 60s.
   - Choose a pair considering off-target hits, amplicon size and exon junctions.
   - Hand-copy `amplicon_start` (forward primer start) and `amplicon_end` (reverse primer start) into an `Amplicons` sheet.
3. **`process_amplicons(xlsx)`** fetches the cDNA by RefSeq ID and cuts `substr(cdna, amplicon_start, amplicon_end)`. It writes:
   - The Excel file.
   - FASTA, with names `GENE_chrN_pos`.
   - Amplicon BED: 10 columns, `0..len`, strand `+`.
   - Site BED: 10 columns, `rel, rel+1`, where `rel = tp − amplicon_start` (**0-based**).
   - A warning if the base at the site isn't T.

---

## 2. Decisions (source: Becca, 2026-09-28)

| # | Topic | Decision |
|---|---|---|
| P-D1 | Automation | **Both modes.** *Auto* mode: local Primer3 plus local specificity search, fully offline. *Import* mode: the user supplies chosen primers (e.g., from a manual PrimerBLAST run) and the tool does everything else. |
| P-D2 | Transcript source | **MANE Select** (the Ensembl and RefSeq sequences are identical), from a pinned release. Fall back to Ensembl canonical, flagged, when a gene has no MANE transcript. |
| P-D3 | Output oligos | **Full ready-to-order oligos.** Forward = `TTTCTGTTGGTGCTGATATTGCG` + forward GSP. Reverse = `ACTTGCCTGTCGCTCTATCTTC` + `N`×10 (UMI) + reverse GSP. Also output the GSPs alone. |
| P-D4 | Integration | **Standalone tool** (`nma-design` CLI and Python package) in the NanoModAmp repo. It is not a Nextflow workflow. Its FASTA/BED outputs are direct inputs to the pipeline. |
| P-D5 | Multiplexing | **Optional** panel-level checks (`--multiplex`): cross-dimers, cross-amplification, and similar Tm across the panel. |
| P-D6 | Organisms | **Human GRCh38 + MANE by default.** Other species through a pinned Ensembl release (canonical transcripts) plus a user-supplied genome and annotation. |
| P-D7 | Site scope | **CDS, UTRs and non-coding RNA exons.** Intronic sites are reported and skipped. |
| P-D8 | Exon junctions | **No preference.** Report whether each primer or amplicon spans a junction, and use it only as the last tie-breaker (configurable). |
| P-D9 | Nearby sites | **Merge** sites on the same transcript that are within **300 nt** of each other into one amplicon (§4.3). |
| P-D10 | Other Ψ sites under primers | **No constraint.** Only the ±50 nt window around target sites applies. |
| P-D11 | Legacy parity | Import mode, given the paper's primers, must reproduce the paper's amplicon FASTA and BEDs exactly. The only allowed differences come from documented register entries (§8), such as a MANE vs Ensembl-canonical difference. |

---

## 3. Architecture

```
design/
  pyproject.toml                  # package: nanomodamp_design; entry point: nma-design
  src/nanomodamp_design/
    refdata.py      # download/verify pinned reference bundles; build indexes
    sites.py        # parse sites; T-strand check; transcript mapping (MANE); region class
    cluster.py      # merge sites ≤ 300 nt (P-D9)
    design.py       # primer3-py candidate generation with window constraints
    specificity.py  # off-target search (BLAST+ blastn-short) + in-silico PCR
    multiplex.py    # optional panel checks
    select.py       # ranking + selection
    importer.py     # import mode (manual primers → coordinates)
    outputs.py      # FASTA, BEDs, TSV/XLSX, order sheet, report
    cli.py
  containers/Dockerfile           # python, primer3-py, BLAST+, pyfaidx, pandas, openpyxl
  tests/                          # pytest + fixtures (mini genome/transcriptome)
docs/design.md                    # user guide + methods
```

**CLI:**

```
nma-design refdata  --species human --release <MANE vX / Ensembl N> --outdir refs/
nma-design annotate --sites sites.txt --refdata refs/ --out run/          # step 1 only
nma-design design   --sites sites.txt --refdata refs/ --config design.yaml --out run/   # auto mode
nma-design import   --primers chosen.tsv --refdata refs/ --out run/        # import mode
nma-design primerblast-sheet --run run/     # emit PrimerBLAST-ready inputs (supports manual use)
```

The tool is fully offline once `refdata` has run. There are **no live biomaRt or Ensembl queries at design time** (P-02).

---

## 4. Algorithm specification

### 4.1 Reference data (`refdata`)
- **Human:**
  - GRCh38 primary assembly FASTA.
  - MANE Select GTF + transcript FASTA (RefSeq and Ensembl) + summary table, at a pinned MANE release.
  - The full RefSeq and Ensembl transcriptomes at pinned releases, used **only for specificity**, to catch paralogs and non-MANE isoforms.
- **Other species:** a user-supplied genome FASTA + GTF, or a pinned Ensembl release. Canonical = the Ensembl canonical tag.
- Record the URLs, release numbers and md5s in `refs/manifest.json`. Build `.fai` indexes and BLAST databases.
- The agent must look up the current MANE and Ensembl releases, pin them, and record them in `docs/DECISIONS.md` **[GATE P0]**.

### 4.2 Site annotation (port of `process_coordinates`, with fixes)
Input: a text or TSV file of sites, `chrN:pos` (1-based, GRCh38, UCSC-style names). An optional `name` column and an optional `strand` column may be included.

For each site:

1. **Strand:** look up the genomic base. If it is T, the strand is `+`; if it is A, the strand is `−`. Otherwise the status is `no_T`: report it and skip. If `strand` was given, it must agree, or the status is `strand_mismatch`.
2. **Transcripts:** find transcripts on that strand whose **exons contain the site**. Take the MANE Select transcript if there is one. Otherwise take the Ensembl canonical transcript and set `transcript_source = ensembl_canonical`. If no exon contains the site but it lies within a transcript, the status is `intronic`: report it and skip (P-D7). If several genes qualify (overlapping genes on the same strand), choose by MANE, then canonical, then protein-coding, and list the alternatives in `alt_transcripts`.
3. **Transcript position:** compute the 1-based `transcript_position` from the exon structure of the **same** transcript whose sequence will be used (P-01). For minus-strand genes, order the exons from 5′ to 3′.
4. **Region class:** `5UTR`, `CDS`, `3UTR` or `ncRNA_exon`.
5. **Check:** the transcript sequence at `transcript_position` must be `T`. If not, it is a **hard error** for that site (P-05).
6. **Output columns** (`sites_annotated.tsv`):
   ```
   coordinate, chr, pos, strand, gene_name, refseq_mrna, ensembl_transcript_id, transcript_source,
   transcript_length, position_type, transcript_position, status, alt_transcripts
   ```
   The first ten names are kept from the legacy table for compatibility.

### 4.3 Merge nearby sites (P-D9)
- Group sites by transcript and sort them by `transcript_position`.
- Build clusters greedily. A site joins the current cluster if its distance from the cluster's **first** site is ≤ `merge_distance` (default 300). Otherwise it starts a new cluster. This means each cluster spans at most 300 nt. **[GATE P1: confirm "span ≤ 300" rather than chaining]**
- Amplicon ID: `GENE_chrN_pos`, using the cluster's **first** site (5′-most in the transcript). This keeps the legacy naming, which downstream code and the MPRA naming rely on. Every member site is listed in `amplicon_sites.tsv` and gets its own site-BED row.
- `merge_distance = 0` turns merging off (legacy behavior: one amplicon per site).

### 4.4 Primer design (auto mode)
The template is the transcript sequence. For a cluster covering positions `[s_first, s_last]`:

- **Window constraints** (legacy ±50): the forward primer must lie entirely within `[1, s_first − 50]`, and the reverse primer entirely within `[s_last + 50, len]`. If either side doesn't have room for a primer of minimum length, the status is `insufficient_flank`, and the report includes the available lengths. The legacy docs said to "adjust accordingly"; the tool relaxes the flank down to `min_flank` (default 25) and records that it did.
- **Primer3 settings** (`design.yaml` defaults, all configurable; confirm at **[GATE P1]**):

  | Setting | Default | Source |
  |---|---|---|
  | `PRIMER_PRODUCT_SIZE_RANGE` | 200–1000 | legacy |
  | `PRIMER_MAX_TM` | 72.0 | legacy |
  | `PRIMER_OPT_TM` | 68.0 | legacy ("upper 60s") |
  | `PRIMER_MIN_TM` | 62.0 | proposed |
  | `PRIMER_PAIR_MAX_DIFF_TM` | 3.0 | proposed (PrimerBLAST default) |
  | `PRIMER_MIN/OPT/MAX_SIZE` | 18 / 24 / 30 | proposed (the higher Tm needs longer primers) |
  | GC, poly-X, self-complementarity | Primer3 defaults | proposed |
  | `PRIMER_NUM_RETURN` | 20 | proposed |

- Tm and dimer checks use the **gene-specific part only**. Hairpin and dimer checks are also run on the **full oligos**, with tails and `N10` modelled as a neutral placeholder. Report both.

### 4.5 Specificity (a PrimerBLAST-equivalent)
- Align each GSP with BLAST+ `blastn-short` against (a) the full transcriptome and (b) the genome.
- **Off-target rule** (PrimerBLAST defaults, configurable): a hit counts as a potential off-target binding site unless the primer has ≥ 2 total mismatches to it, including ≥ 2 within the last 5 nt at the 3′ end. Hits with ≥ 6 mismatches are ignored.
- **In-silico PCR:** an off-target *product* exists when a forward and a reverse binding site (either primer, any orientation pairing) fall on the same sequence, facing each other, within `max_offtarget_product` (default 4000 nt).
- Report, per pair:
  - `n_offtarget_products_transcriptome` and `n_offtarget_products_genome`.
  - The list of off-target products.
  - `n_binding_sites` for each primer.
  - Isoform hits: products on other isoforms of the *same gene* with the same insert sequence are marked `same_gene_isoform`. They are not penalized.

### 4.6 Ranking and selection
Candidates are sorted lexicographically, from most to least important:

1. Fewest off-target products (transcriptome, then genome).
2. Primer3 pair penalty.
3. Product size closest to `preferred_product_size` (default: the product-range midpoint; **[GATE P1]** Becca may give a preferred size).
4. Junction flag, as the last tie-breaker (P-D8, off by default).

Choose the top pair. Keep the top `n_report` pairs (default 5) in `candidates.tsv` for manual review.

If no pair has zero off-target products, still select the best pair, but set its status to `offtarget_warning`.

### 4.7 Optional multiplex checks (`--multiplex`, P-D5)
Using the selected full oligos across the whole panel:

- **Cross-dimers:** Primer3 heterodimer ΔG for every pair of primers, flagged below a threshold (default −9 kcal/mol).
- **Cross-amplification:** in-silico PCR with every primer combination against the transcriptome and genome.
- **Tm spread:** the panel range of GSP Tm, with a warning if it exceeds 5 °C.
- **Optional iterative re-selection:** swap in the next-best candidates for flagged amplicons, up to `max_iterations`.

Output: `multiplex_report.tsv`, plus a suggested split into compatible pools when conflicts remain.

### 4.8 Import mode (P-D1)
- Input is `chosen.tsv` with the columns `amplicon_id` or `coordinate`, `fwd_gsp`, `rev_gsp`. The legacy `Amplicons` Excel sheet is also accepted, and read for its `refseq_mrna`, `amplicon_start` and `amplicon_end` columns.
- **Locate** each primer on the MANE transcript by exact match: forward on the sense strand, reverse as its reverse complement. The amplicon runs from the start of the forward primer to the end of the reverse primer's binding site, inclusive.
  - If a primer is not found, or found more than once, that is an error.
  - If only legacy coordinates are given, cut the sequence from them directly, and verify the T at every site.
- Then run §4.5, the optional §4.7 and §4.9. Specificity is reported only; nothing is selected.

### 4.9 Outputs (both modes)
| File | Content |
|---|---|
| `sites_annotated.tsv` | §4.2 |
| `amplicon_sites.tsv` | amplicon_id ↔ member sites, with the relative (0-based) position of each |
| `primers.tsv` | Per amplicon: fwd/rev GSP, positions, Tm, GC, penalty, product size, junction flags, off-target counts, status |
| `candidates.tsv` | Top-N pairs per amplicon (auto mode) |
| `order_oligos.tsv` | `name, sequence, scale, purification, notes`. The forward is `<id>_F` = tail + GSP; the reverse is `<id>_R` = tail + `NNNNNNNNNN` + GSP. The notes say to order `N` as a machine-mixed or hand-mixed random base (vendor-specific). Vendor columns are configurable; the default is an IDT-compatible layout. |
| `amplicons.fasta` | `>GENE_chrN_pos` + amplicon sequence (**the reference FASTA for the endogenous pipeline**) |
| `amplicons.bed` | 10 columns, as in the legacy file: `id 0 len GENE 0 + NA gene NA NA` |
| `amplicons_position.bed` | 10 columns, one row per site: `id rel rel+1 GENE 0 + NA gene NA NA` (standard 0-based BED) |
| `amplicons.xlsx` | The legacy-compatible workbook, with an `Amplicons` sheet plus all columns above |
| `primerblast_inputs.tsv` | Per site/cluster: RefSeq ID, "forward primer to", "reverse primer from", product range and Tm settings (supports manual PrimerBLAST) |
| `design_report.html` | Summary: status counts, per-amplicon diagrams (site positions, primers, junctions), off-target tables, multiplex results |
| `run_manifest.json` | Tool version, config, refdata manifest, timestamps |

**Pipeline hand-off:** `amplicons.fasta` → `--fasta`, and `amplicons.bed` or `amplicons_position.bed` → `--bed`, with `bed_coordinates = bed0`. See §6 about R-15.

---

## 5. Tests

**Fixtures:** a mini "genome" of three 5-kb contigs and matching GTF/transcript FASTA. It contains:

- A plus-strand gene and a minus-strand gene, each with ≥ 3 exons.
- A non-coding RNA.
- A gene with two isoforms (MANE plus an alternative).
- A paralog with 1 and 3 mismatches in the primer regions.
- Two overlapping genes on the same strand.

Also build a small real-data fixture from 3–5 paper sites, using the pinned human MANE extracted to a small subset (check the size).

**Unit tests (pytest):**

1. `chr:pos` parsing; A→`−`, T→`+`, and G/C → `no_T`.
2. Transcript position is correct at exon boundaries, for both strands. Include the first and last base of an exon.
3. UTR, CDS and ncRNA classification.
4. An intronic site is reported and skipped.
5. MANE is preferred, with fallback to canonical, and the flag is set.
6. The T check fails loudly when the sequence and coordinates disagree.
7. **Merging:** sites 0, 150 and 300 nt apart form one cluster; sites 301 apart form two; `merge_distance=0` gives one amplicon per site; the cluster ID uses the first site.
8. **Windows:** no primer overlaps `[s_first−49, s_last+49]`. An insufficient flank gives the right status, and the relaxation to `min_flank` is recorded.
9. **Specificity:** the 1-mismatch paralog is detected as an off-target product; the 3-mismatch-in-3′-end paralog is ignored per the rule; same-gene isoform products are labelled and not penalized.
10. The ranking order is deterministic.
11. Full oligos are assembled exactly as `tail + [N10] + GSP`. The GSP Tm excludes the tails.
12. **Multiplex:** a planted heterodimer is flagged, and a planted cross-amplicon is detected.
13. **Import mode:**
    - Primers are located correctly on both strands.
    - An amplicon from legacy coordinates is extracted identically.
    - Not-found and multi-found primers raise errors.
14. **Outputs:**
    - BED rows are 0-based.
    - `amplicons_position.bed` start + 1 indexes a `T` in `amplicons.fasta` for every row. This is the key cross-file invariant.
    - FASTA names match the BED chrom names.
    - Outputs pass the NanoModAmp input schema validation.
15. Identical inputs and config give byte-identical outputs.

**Legacy parity (golden) test, P-D11 [GATE P2]:**

- Becca provides the paper's endogenous design workbook (the `Amplicons` sheet, including primers if available) and the paper's amplicon FASTA and BEDs.
- `nma-design import` must reproduce the paper's FASTA sequences and both BEDs exactly. Every difference must be explained in a register entry (§8), for example when the paper's Ensembl-canonical cDNA differs from MANE in UTR length (P-01).
- Separately, run auto mode on the same sites and report, per site, whether it found a valid pair, how it compares with the paper pair (overlap and product size), and the off-target counts. This part is for information only.
- Store the fixture in `design/tests/data/paper/` (small text files only).

**Integration:** run `nma-design design` on the fixture, then run the NanoModAmp synthetic simulator (WP1) on the resulting FASTA/BED, then run the pipeline `-profile test`. The pipeline must accept the files and report the designed sites.

---

## 6. Links to the other plans

- **R-15 (BED convention, `plan.md` §9):** the legacy design code writes **standard 0-based** site BEDs (`rel, rel+1`, where `rel` = 0-based offset). The amplicon BED starts at 0. `bam_counts_fast.R` treats the BED start as 1-based (its `+1` is commented out). If the paper's endogenous BEDs came from this tool, the missing `+1` in counting is a bug: for a single-site BED, it pulls in the upstream base. WP5 should use this as the leading hypothesis and confirm it with the golden data (golden plan §4.9).
- **`sample_name.R` (D26):** the original exists at `rmrodell/BIDamplicon/legacy/sample_name.R`. `BIDamplicon/pipeline/BIDdetect.sh` and `bam_counts_fast.R` are identical to the paper versions. The golden-package agent should use this file instead of reconstructing it.
- **Endogenous reference (D6):** in practice, the reference for the endogenous pipeline is `amplicons.fasta`, i.e. transcript-derived amplicon sequences.

---

## 7. Work breakdown (WP9)

| Step | Task | Depends on | Done when |
|---|---|---|---|
| 9.0 | Package scaffold, CLI skeleton, container, CI job (pytest) | — | CI green |
| 9.1 | `refdata`: pin MANE/Ensembl/GRCh38, download, verify, index | 9.0 | **[GATE P0]** releases approved |
| 9.2 | Fixtures (mini genome) | 9.0 | Fixtures committed |
| 9.3 | Site annotation (§4.2) + merging (§4.3) | 9.1, 9.2 | Tests 1–7 pass |
| 9.4 | Primer3 design (§4.4) | 9.3 | Tests 8, 11 pass; **[GATE P1]** defaults approved |
| 9.5 | Specificity + in-silico PCR (§4.5), ranking (§4.6) | 9.4 | Tests 9–10 pass |
| 9.6 | Import mode (§4.8), outputs (§4.9), report | 9.3 | Tests 13–15 pass |
| 9.7 | Multiplex (§4.7) | 9.5 | Test 12 passes |
| 9.8 | Legacy parity + integration | 9.6, WP1 | **[GATE P2]** parity report approved |
| 9.9 | Docs (`docs/design.md`): quick start, both modes, PrimerBLAST hand-off, parameters, methods, limitations | 9.6 | Reviewed |

---

## 8. Change register for the design tool (add to `CHANGE_REGISTER.md` with the prefix `P-`)

| ID | Legacy behavior | Resolution | Status |
|---|---|---|---|
| P-01 | Position computed from Ensembl exons; primers designed on the RefSeq mRNA; amplicon cut from Ensembl cDNA. Offsets arise when UTRs differ | One MANE transcript for everything | APPROVED |
| P-02 | Live, unpinned Ensembl (biomaRt) | Offline, pinned reference bundle | APPROVED |
| P-03 | `canonical_transcripts[1,]` and `cdna[1]` pick arbitrarily among multiple matches | Deterministic MANE → canonical → biotype order, with alternatives reported | APPROVED |
| P-04 | Rows with `no T` have fewer columns, so `rbind` fails inside `tryCatch` and sites are **dropped silently** (after that, every later site may fail too) | A fixed schema with a `status` column for every site | APPROVED |
| P-05 | A non-T at the site only produces a warning | Hard error per site | APPROVED |
| P-06 | Errors swallowed by `tryCatch` warnings | Explicit per-site `status` values and a summary | APPROVED |
| P-07 | Docs say "CDS only", but the code accepts any exon | Explicit scope: CDS, UTR, ncRNA (P-D7) | APPROVED |
| P-08 | Manual PrimerBLAST, Excel hand-entry | Auto mode plus import mode (P-D1) | APPROVED |
| P-09 | One amplicon per site | Merge ≤ 300 nt (P-D9); `merge_distance=0` gives legacy behavior | APPROVED |
| P-10 | Excel is the primary interchange format | TSV is primary; XLSX is kept for compatibility | APPROVED |
| P-11 | Human only | Human MANE by default, other species configurable (P-D6) | APPROVED |

---

## 9. Gates and open items

| Gate | Becca decides / provides |
|---|---|
| P0 | Pinned MANE, Ensembl and GRCh38 releases |
| P1 | Primer3 defaults in §4.4 (especially min Tm, primer lengths, preferred product size); merge rule "cluster span ≤ 300 nt"; the vendor order-sheet format; whether the reverse-primer N10 should be ordered as machine-mixed |
| P2 | Parity report against the paper's design files; provides the paper workbook, FASTA and BEDs |

Open questions to confirm at P1:

- Is the reverse primer also the RT primer (i.e., gene-specific RT with UMI), or is RT done separately? This affects whether dimer checks should include RT conditions.
- Should primers avoid common SNPs (e.g., dbSNP common variants in the primer 3′ end)? This is off by default.
- Is there a maximum panel size for multiplex mode?
