# Sample sheet (plan §5.1)

CSV with a header row. Schema: `assets/schema_input.json` (validated with nf-schema), plus
checks in `readSamplesheet()` (`subworkflows/local/utils_nfcore_nanomodamp_pipeline`).

| Column | Required | Rule |
|---|---|---|
| `sample_id` | yes | Unique; matches `[A-Za-z0-9._-]+`. |
| `fastq` | yes | A `.fastq.gz` file, or a directory whose `*.fastq.gz` files are concatenated in lexical order. Relative paths are resolved against the sample sheet's directory. |
| `treat` | yes | Exactly `input` or `BS` (D14). |
| `rep` | yes | Integer or string; batch-paired replicate (D17). |
| others | no | Metadata, carried in sample-sheet order into `counts_merged.tsv` and usable in `subset` and `factor`. |

Errors: duplicate `sample_id`; missing file; empty FASTQ directory; `treat` not `input`/`BS`;
a `fastq` that is neither `.fastq.gz` nor a directory. Blank lines are ignored with a warning.

Example: [examples/samplesheet.csv](examples/samplesheet.csv).
