#!/bin/bash
# Source after REF_ENDO / REF_POOL1 are set. Copies both references into the build dir with a
# freshly built .fai and repoints REF_ENDO / REF_POOL1 at the copies (inputs are never modified).
# Why: pool1_cleaned_noadapters.fasta.fai (2025-11-07) predates a CRLF->LF rewrite of the FASTA
# (2026-01-26), so its offsets are wrong and Rsamtools would read wrong ref bases/kmers.
STAGED_REFS=${GOLDEN_DIR:?}/refs
mkdir -p "$STAGED_REFS"
module load biology samtools/1.16.1 >/dev/null 2>&1
for v in REF_ENDO REF_POOL1; do
  src=${!v}; dst=${STAGED_REFS}/$(basename "$src")
  cmp -s "$src" "$dst" || cp "$src" "$dst"
  samtools faidx "$dst"
  printf -v "$v" '%s' "$dst"
done
export REF_ENDO REF_POOL1
