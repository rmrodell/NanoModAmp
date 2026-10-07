# Counts (plan §5.3)

**Per sample:** `counts/per_sample/<sample>.counts.tsv`, tab-separated. Columns, exactly and in
this order (paper-compatible):

```
chr pos gene totalReads A.count C.count G.count T.count Deletion.count Insertion.count ref kmer strand delrate
```

Types: `pos` and all count columns are integers; `delrate` = `Deletion.count / totalReads`;
`totalReads = A+C+G+T+Deletion+Insertion` (D9); only rows with `totalReads > min_coverage` (D10).
`kmer` is NA within 2 nt of either end. `ref` is `T` only unless `count_all_bases` (D12).

**Merged:** `counts/counts_merged.tsv` = `sample_id` + all sample-sheet metadata columns (in
sample-sheet order, excluding `fastq`) + the per-sample columns above. Tables are joined **by
column name, never by position** (R-06).

**Merging runs (D32):** with `--input_counts`, one or more `counts_merged.tsv` tables from earlier
runs are merged into one table with the same layout, written as `counts/counts_merged.tsv`, plus
`counts/merge_sources.tsv` (`sample_id`, `source_table`). Rules: columns are matched by name (R-06);
metadata columns are the union, in first-seen order (missing values become NA, with a warning);
the count and site columns must be identical across tables; a `sample_id` may appear in only one
table; sites missing in a sample are simply absent (no fill). Example:
[examples/merge_sources.tsv](examples/merge_sources.tsv).

**Failed regions:** `counts/per_sample/<sample>.failed_regions.tsv` with `chr, start, end, gene, error`
(always written; header only when nothing failed; R-12).

Examples: [examples/sample.counts.tsv](examples/sample.counts.tsv), [examples/counts_merged.tsv](examples/counts_merged.tsv).
