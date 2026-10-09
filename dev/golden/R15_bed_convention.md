# R-15: BED coordinate convention used for the paper counts

**Status: decided by Becca at G1-e, 2026-10-06:** bed0 default; golden targets.bed converted to 0-based; 20250418 window-start rows not intended (paper artifact). Written 2026-10-06 (golden plan §4.9). All inputs were read-only. Reference bases were read from freshly indexed copies in `/scratch/users/rodell/NanoModAmp_build/golden/refs/r15/`, because the on-disk `pool1_cleaned_noadapters.fasta.fai` is stale.

## 1. How the legacy counters read a BED

All three counters behave the same way here: `dev/legacy/endogenous/20250418/bam_counts.R`, `.../20251022/bam_counts.R` and `mpra_*/bam_counts_fast.R`.

```r
bedFile <- bedFile[, c(1,2,3,4,6)]; colnames = chr, start, end, gene, strand
# bam_counts_fast.R only:
# # BED start is 0-based, GenomicRanges is 1-based. Add 1 to start.
# bedFile[, start := start + 1]          <- commented out
GRanges(seqnames = chr, ranges = IRanges(start = start, end = end), strand = strand)
```

- BED `start` and `end` go straight into a 1-based, closed `IRanges`. The pileup is restricted to that range, so the legacy code reports **every 1-based position from `start` to `end` inclusive**.
- A standard 0-based, half-open BED would instead mean positions `start+1` to `end`.
- The only difference between the two readings is therefore one position: `start` itself. Legacy counts it; standard BED excludes it.
- The `ref == "T"` filter then drops that extra position unless the base there is T.
- Checks confirming the counters only report positions inside the range:
  - both shipped MPRA delpos tables (4,574 and 42,178 rows) have `pos` equal to the BED `end` in every row;
  - all endogenous table positions fall inside their BED rows.

## 2. Evidence per BED

"T@start" means the 1-based base at `start` is T; "T@start+1" means the base at `start+1` is T. A designed Ψ site is always a T.

### Single-site BEDs

| BED (used for) | Rows | start = end | T@start | T@start+1 | Convention | Paper positions vs BED |
|---|---|---|---|---|---|---|
| `pool1_cleaned_delpos_noadapters.bed` (MPRA in vitro and in cellulo, delpos) | 775 | 775 | **775** | 277 | 1-based site, `start = end = pos` | `pos = end` in all 46,752 rows; no extra position |
| `set3_delpos.bed` (20251022) | 7 | 6 | 6 | 1 | **Mixed**: 6 rows 1-based (`start = end`); STIM1 `260 261` is standard 0-based | STIM1: only 261 reported. 260 is G, so the extra position was filtered |
| `endoG1G2_delpos.bed` (golden, F2) | 10 | 6 | 6 | 5 | **Mixed**: 6 rows 1-based (`start = end`); HP1BP3 `324 325`, INTS10 `311 312`, SMOC1 `166 167`, ERCC8 `281 284` are standard 0-based | Only the start=end rows are golden targets (RHBDD2 287, HDAC6 392, STIM1 261) |
| `4psi.bed` | 4 | 4 | 4 | 0 | 1-based site | — |

The 1-based reading is confirmed by the bases:

| Row | 1-based base at `start` | base at `start+1` |
|---|---|---|
| RHBDD2 287 287 | T | A |
| HDAC6 392 392 | T | A |
| STIM1 261 261 | T | A |
| 4psi 35/70/104/136 | T | A / G / G / C |

The 0-based rows also put the T at `end`:

| Row | base at `start` | base at `start+1` = `end` |
|---|---|---|
| HP1BP3 324 325 | C | T |
| STIM1 260 261 (set3) | G | T |

Under a standard reading, `start = end` is a zero-length interval and would count nothing.

### Window BEDs (several positions per row)

| BED (used for) | Example row | Bases at start / start+1 / end | Legacy positions | Standard positions | Paper table |
|---|---|---|---|---|---|
| `20250404_endoBIDamplicon.bed` (**20250418 paper counts**) | RHBDD2 `285 287`, HDAC6 `390 392`, ACTL6A `87 89`, EIF5 `145 147` | T / G / T in all four | 285–287 | 286–287 | Reports **both 285 and 287**. 78 of the 248 rows are at the window `start` |
| `20250508_293Thigh_PUS7dep_amplicons_position.bed` (20250522, 20250606) | UBFD1 `141 147`, IGF2BP1 `475 481`, CHP1 `400 406` | `start` is non-T in 8 of 9 rows (EIF5 145 is T but EIF5 isn't in those tables) | p−4 … p+2 | p−3 … p+2 | No row at `start`; table positions (e.g. UBFD1 142–145, CHP1 401–404) fall inside both readings |

### Full-amplicon BEDs

| BED | start values | `end` vs sequence length | Effect |
|---|---|---|---|
| `pool1_cleaned_full_noadapters.bed` | all `1` | `130` = length | Legacy counts position 1; a standard reading drops it. The paper full tables have pos = 1 rows: in vitro 1,359 rows (228 oligos), in cellulo 12,174 rows (229 oligos) |
| `set3_fulllength.bed`, `20250404_endoBIDamplicon_full.bed` | `1` | STIM1 `520` for a 521 nt sequence (last base never counted); others = length | Position 1 is counted only if T (ACTL6A and EIF5: no) |
| `endoG1G2_full.bed` | `0` for 7 rows, `1` for HDAC6, RHBDD2, STIM1 | STIM1 `520` of 521 | `IRanges(start = 0)` covers a nonexistent position 0, so the effect is the same as starting at 1 |

Side note for WP3, outside R-15: at pos 1 the paper full tables write `kmer = ""` rather than NA.

## 3. Conclusions

1. **The paper's single-site BEDs are 1-based**, with `start = end = pos`: pool1 delpos, 4psi, and the `start = end` rows of `endoG1G2_delpos` and `set3_delpos`. The legacy code (no +1) counts exactly the designed T for them. All MPRA delpos counts and all golden-target counts are therefore at the correct position.
2. **Some rows in `endoG1G2_delpos` and `set3_delpos` are standard 0-based**: HP1BP3, INTS10, SMOC1, ERCC8, and STIM1 in set3. Legacy also counts the upstream base for these, but it is non-T in every case, so the paper tables have no spurious rows from them. The files mix conventions.
3. **The paper counts contain one possible spurious-position case: 20250418.** The window BED `20250404_endoBIDamplicon.bed` (`285 287` etc.) gave rows at the window `start` (285, 390, 87, 145): 78 of 248 rows. Under a standard reading those positions fall outside the windows. Under a 1-based inclusive reading they are intended.
   - The designed sites are at `end` (RHBDD2 287, HDAC6 392, matching `endoG1G2_delpos`).
   - The base at `start` is T in all four windows, so this could be deliberate.
   - **Becca needs to say which was intended.** It doesn't affect the golden package: V1 compares only 287 and 392.
4. **The paper full-amplicon MPRA tables include oligo position 1** (start = 1), which a standard 0-based reading would exclude. That is about 0.8% of rows, at the oligo edge, and none are delpos sites.
5. Standard 0-based BEDs, as `nma-design` will write (`primer_plan.md`: site BED `rel, rel+1`), read correctly only if the pipeline adds +1. Under the legacy code, `start = end` BEDs read correctly and 0-based BEDs pick up one extra upstream position, which is counted only if it is a T.

## 4. Recommendation for `bed_coordinates` (plan §3)

- **Default `bed0`** (standard: positions `start+1 … end`). This matches every BED tool and the planned `nma-design` output, and gives correct results for all 0-based rows above.
- **`one_based_start`** (legacy: positions `start … end`) for the paper's existing BEDs. It reproduces the paper positions exactly, including the 20250418 window starts and MPRA position 1.
- **Validation:**
  - Under `bed0`, any row with `start = end` (a zero-length interval) should be an **error**, with a message saying the file looks like a 1-based single-site BED and pointing to `bed_coordinates = one_based_start`. This catches the paper-style delpos files instead of silently counting nothing.
  - Do not add an `auto` mode, because the paper files mix conventions.
- **Golden package (`targets.bed`)**, two options:
  - **A (recommended):** keep the paper rows unchanged (golden plan §4.7) and declare `bed_coordinates = one_based_start` in `conf/test_golden.config` and `analyses`/README. All selected rows are `start = end` (endogenous 287, 392, 261; MPRA delpos 65-style sites), so the counts are identical to the paper's.
  - **B:** convert `targets.bed` to standard 0-based (`pos−1, pos`) and use the default `bed0`. Counts are again identical, but the file no longer matches the paper bytes. Record the conversion in provenance.

**Becca decides at G1-e:**
1. the default (recommended `bed0`);
2. option A or B for the golden `targets.bed`;
3. whether the 20250418 window-start positions (285, 390, 87, 145) were intended (relevant only to WP7 concordance against the full 20250418 paper table, not to the golden V1).
