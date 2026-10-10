#!/bin/bash
# Golden plan §4.2 check L1: rerun BAM read IDs vs paper BAM read IDs, per sample and
# selected chr. Compares the set of read names (L1 criterion) and, for information,
# the set of alignments (QNAME, FLAG, RNAME, POS, CIGAR) and the multiset of primary
# (RNAME, POS, strand) positions. umi_tools dedup picks a random read among ties (no seed was
# set in the paper runs), so dedup BAMs can keep the same groups with different representatives:
# positions_equal = yes with ids_equal = no means "same dedup groups, different representatives".
# Both the deduplicated BAM and the pre-dedup BAM are checked where both exist.
#
# Usage: l1_compare.sh   (env from submit_l1.sh: L1_DIR, CHR_ENDO, CHR_MPRA, PAPER_*)
# Output: $L1_DIR/L1_results.tsv
set -euo pipefail
: "${L1_DIR:?}"
module load biology samtools/1.16.1 >/dev/null 2>&1

RES="${L1_DIR}/L1_results.tsv"
TMP=$(mktemp -d "${L_SCRATCH:-/tmp}/l1cmp.XXXX")
echo -e "experiment\tsample\tbam\tchr\tpaper_ids\trerun_ids\tshared_ids\tpaper_only\trerun_only\tids_equal\talignments_equal\tpositions_equal" > "$RES"

# One streaming pass per BAM (no index needed, nothing written next to input BAMs):
# writes $TMP/<tag>.<chr>.aln = sorted "QNAME FLAG RNAME POS CIGAR" for each selected chr.
split_bam() {  # bam tag chr_file
  rm -f "$TMP/$2".*
  samtools view "$1" | awk -F'\t' -v OFS='\t' -v d="$TMP" -v t="$2" -v cf="$3" \
    'BEGIN{while((getline c < cf)>0) if(c!="") keep[c]=1}
     ($3 in keep){print $1,$2,$3,$4,$6 > (d"/"t"."$3".aln")}'
  while read -r c; do [ -n "$c" ] || continue
    touch "$TMP/$2.$c.aln"; sort -o "$TMP/$2.$c.aln" "$TMP/$2.$c.aln"
    cut -f1 "$TMP/$2.$c.aln" | sort -u > "$TMP/$2.$c.ids"
    awk -F'\t' '!and($2,2308){print $3, $4, and($2,16)}' "$TMP/$2.$c.aln" | sort > "$TMP/$2.$c.pos"
  done < "$3"
}

compare() {  # exp sample label paper_bam rerun_bam chr_file
  local exp=$1 s=$2 lab=$3 pb=$4 rb=$5 chrs=$6 c f
  for f in "$pb" "$rb"; do
    if [ ! -s "$f" ]; then echo -e "${exp}\t${s}\t${lab}\tALL\tNA\tNA\tNA\tNA\tNA\tMISSING:$(basename "$f")\tNA\tNA" >> "$RES"; return; fi
  done
  split_bam "$pb" p "$chrs"; split_bam "$rb" r "$chrs"
  while read -r c; do
    [ -n "$c" ] || continue
    local np nr ns ae=no
    np=$(wc -l < "$TMP/p.$c.ids"); nr=$(wc -l < "$TMP/r.$c.ids"); ns=$(comm -12 "$TMP/p.$c.ids" "$TMP/r.$c.ids" | wc -l)
    cmp -s "$TMP/p.$c.aln" "$TMP/r.$c.aln" && ae=yes
    local pe=no; cmp -s "$TMP/p.$c.pos" "$TMP/r.$c.pos" && pe=yes
    echo -e "${exp}\t${s}\t${lab}\t${c}\t${np}\t${nr}\t${ns}\t$((np-ns))\t$((nr-ns))\t$([ "$np" = "$ns" ] && [ "$nr" = "$ns" ] && echo yes || echo no)\t${ae}\t${pe}" >> "$RES"
  done < "$chrs"
}

# endogenous 20250418 (paper counts came from the pre-dedup sorted/ BAMs)
while read -r fq; do s=$(basename "$fq" .fq)
  compare endo_20250418 "$s" pre_dedup "${PAPER_0418}/sorted/${s}_UMI_sort.bam" "${L1_DIR}/endo_20250418/${s}/minimap2/sorted/${s}_UMI_sort.bam" "$CHR_ENDO_0418"
  compare endo_20250418 "$s" dedup "${PAPER_0418}/dedup/${s}_UMI_dedup.bam" "${L1_DIR}/endo_20250418/${s}/minimap2/dedup/${s}_UMI_dedup.bam" "$CHR_ENDO_0418"
done < "${L1_DIR}/inputs/endo_20250418/fastqs.txt"

# endogenous 20251022 and MPRA in vitro (same script family, same layout)
for exp in endo_20251022 mpra_invitro; do
  if [ "$exp" = endo_20251022 ]; then pd=$PAPER_1022; chrs=$CHR_ENDO_1022; else pd=$PAPER_INVITRO; chrs=$CHR_MPRA; fi
  while IFS=: read -r s _; do s=$(echo "$s" | xargs)
    compare "$exp" "$s" pre_dedup "${pd}/${s}/tmp/${s}_mapped_sorted.bam" "${L1_DIR}/${exp}/retained/${s}_mapped_sorted.bam" "$chrs"
    compare "$exp" "$s" dedup "${pd}/deduplicated_bam/${s}.bam" "${L1_DIR}/${exp}/deduplicated_bam/${s}.bam" "$chrs"
  done < "${L1_DIR}/inputs/${exp}/sample_map.txt"
done

# MPRA in cellulo (no UMI, no dedup: final_bam is the mapped, sorted BAM)
for fq in "${L1_DIR}"/inputs/mpra_incell/fq/*.fq; do s=$(basename "$fq" .fq)
  compare mpra_incell "$s" final_no_dedup "${PAPER_INCELL}/${s}.bam" "${L1_DIR}/mpra_incell/final_bam/${s}.bam" "$CHR_MPRA"
done

rm -rf "$TMP"
echo "L1 summary (rows with ids_equal):"
echo -e "experiment\tbam\tids_equal\talignments_equal\tpositions_equal"
awk -F'\t' 'NR>1{k=$1"\t"$3; n[k]++; if($10=="yes") a[k]++; if($11=="yes") b[k]++; if($12=="yes") c[k]++}
  END{for(k in n) printf "%s\t%d/%d\t%d/%d\t%d/%d\n", k, a[k], n[k], b[k], n[k], c[k], n[k]}' "$RES" | sort
