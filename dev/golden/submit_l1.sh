#!/bin/bash
# Golden plan §4.2: prepare inputs for the approved minimal sample sets (G1-b, 2026-10-06)
# and submit the legacy preprocessing reruns plus the L1 comparison.
#
# Usage: submit_l1.sh [--prepare-only]
# Writes only under $L1_DIR (default /scratch/users/rodell/NanoModAmp_build/golden/l1).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export LEGACY_DIR="$(cd "${HERE}/../legacy" && pwd)"
export L1_DIR="${L1_DIR:-/scratch/users/rodell/NanoModAmp_build/golden/l1}"

# Inputs (dev/golden/golden_inputs.yaml; decisions in docs/DECISIONS.md)
BID=/oak/stanford/groups/nicolemm/rodell/BIDamplicon
export REF_ENDO=/home/groups/nicolemm/rodell/fastas/endoG1G2.fasta              # F2
export REF_POOL1=/home/groups/nicolemm/rodell/pool1/pool1_cleaned_noadapters.fasta
FQ_0418=${BID}/endo_CellTypeSpec/20250418_endoBID/20250503_70/fastq
export FQ_1022=${BID}/endo_CellTypeSpec/20251022_endoBID/fastq
export FQ_INVITRO=${BID}/Pool1/invitro/maydata/00_concatenated_fastqs
FQ_INCELL=${BID}/20241114_pool1/fastq
MAP_1022=${BID}/endo_CellTypeSpec/20251022_endoBID/barcodes.txt
MAP_INVITRO=${BID}/Pool1/invitro/barcodes_may.txt
export PAPER_0418=${BID}/endo_CellTypeSpec/20250418_endoBID/20250503_70/deduplication/minimap2
export PAPER_1022=${BID}/endo_CellTypeSpec/20251022_endoBID/mapping
export PAPER_INVITRO=${BID}/Pool1/invitro/maydata
export PAPER_INCELL=${BID}/Pool1/incell/prep/final_bam
CHR_ENDO=/scratch/users/rodell/NanoModAmp_build/endoBID_data/endo_chr.txt
export CHR_MPRA=/scratch/users/rodell/NanoModAmp_build/pool1_chr.txt

# --- inputs: approved sample sets -------------------------------------------
IN=${L1_DIR}/inputs
mkdir -p ${IN}/{endo_20250418,endo_20251022,mpra_invitro,mpra_incell/fq}
# chr lists per run (endo_chr.txt has no trailing newline; normalise)
export CHR_ENDO_0418=${IN}/endo_20250418/chrs.txt CHR_ENDO_1022=${IN}/endo_20251022/chrs.txt
grep -E '^(RHBDD2|HDAC6)_' <(cat "$CHR_ENDO"; echo) > "$CHR_ENDO_0418"
grep -E '^STIM1_' <(cat "$CHR_ENDO"; echo) > "$CHR_ENDO_1022"
# 20250418: reps 1 and 3 (rep 2 unusable, see endogenous_preflight_G1b.md)
ls ${FQ_0418}/*.fq | awk -F/ '{split($NF,a,"_"); if(a[3]=="1"||a[3]=="3") print}' > ${IN}/endo_20250418/fastqs.txt
# 20251022: reps 1 and 2
awk -F: '{split($1,a,"_"); if(a[3]=="1"||a[3]=="2") print}' "$MAP_1022" > ${IN}/endo_20251022/sample_map.txt
# in vitro: May run, all 4 samples
grep -v '^\s*$' "$MAP_INVITRO" > ${IN}/mpra_invitro/sample_map.txt
# in cellulo: P101/P102/P3/P4 x input/BS; HepG2 reps 1,2; 293T reps 1,3
rm -f ${IN}/mpra_incell/fq/*.fq
for f in ${FQ_INCELL}/*.fq; do
  IFS=_ read -r ct vec rep trt <<< "$(basename "$f" .fq)"
  case "$vec" in P101|P102|P3|P4) ;; *) continue ;; esac
  if { [ "$ct" = HepG2 ] && [[ $rep == [12] ]]; } || { [ "$ct" = 293T ] && [[ $rep == [13] ]]; }; then
    ln -s "$f" ${IN}/mpra_incell/fq/
  fi
done

declare -A N=(
  [endo_20250418]=$(wc -l < ${IN}/endo_20250418/fastqs.txt)
  [endo_20251022]=$(wc -l < ${IN}/endo_20251022/sample_map.txt)
  [mpra_invitro]=$(wc -l < ${IN}/mpra_invitro/sample_map.txt)
  [mpra_incell]=$(ls ${IN}/mpra_incell/fq/*.fq | wc -l)
)
for e in "${!N[@]}"; do echo "$e: ${N[$e]} samples"; done
[ "${N[endo_20250418]}" = 16 ] && [ "${N[endo_20251022]}" = 16 ] && [ "${N[mpra_invitro]}" = 4 ] && [ "${N[mpra_incell]}" = 32 ] \
  || { echo "unexpected sample counts (want 16/16/4/32)" >&2; exit 1; }

[ "${1:-}" = "--prepare-only" ] && exit 0

# --- submit ------------------------------------------------------------------
# 16 CPUs / 32G per task as in the paper submitters (the 20250418 run used 15 threads).
mkdir -p ${L1_DIR}/slurm_logs
ids=()
for e in endo_20250418 endo_20251022 mpra_invitro mpra_incell; do
  t=02:00:00; [ "$e" = mpra_incell ] && t=04:00:00
  id=$(sbatch --parsable -p normal --job-name="golden_L1_${e}" --array=1-${N[$e]} \
        --cpus-per-task=16 --mem=32G --time=$t \
        --output="${L1_DIR}/slurm_logs/${e}_%A_%a.out" --error="${L1_DIR}/slurm_logs/${e}_%A_%a.err" \
        --export=ALL "${HERE}/l1_task.sh" "$e")
  echo "$e -> job $id"; ids+=("$id")
done
dep=$(IFS=:; echo "${ids[*]}")
cid=$(sbatch --parsable -p normal --job-name=golden_L1_compare --dependency=afterok:${dep} \
      --cpus-per-task=2 --mem=8G --time=02:00:00 \
      --output="${L1_DIR}/slurm_logs/compare_%j.out" --error="${L1_DIR}/slurm_logs/compare_%j.err" \
      --export=ALL "${HERE}/l1_compare.sh")
echo "compare -> job $cid (afterok:${dep})"
echo "${ids[*]} $cid" > ${L1_DIR}/job_ids.txt
