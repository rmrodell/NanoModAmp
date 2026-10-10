Stale-index fixture (R-31). `ref.fa.fai` was built from a CRLF version of `ref.fa`, which was then
rewritten with LF line endings: names and lengths still match, offsets do not (the pool1 case).
`expected_fresh.fai` is `samtools faidx` of the current `ref.fa`. `ref.fa.gz` is gzip of `ref.fa`.
