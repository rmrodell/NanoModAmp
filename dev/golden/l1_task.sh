#!/bin/bash
# One SLURM array task of golden plan §4.2: rerun the legacy preprocessing for one sample.
# Runs the dev/legacy script that produced the paper data (F1), unchanged, on the inputs
# prepared by submit_l1.sh. Never writes outside $L1_DIR.
#
# Usage (from submit_l1.sh): sbatch --array=1-N l1_task.sh <experiment>
#   experiment = endo_20250418 | endo_20251022 | mpra_invitro | mpra_incell
set -eo pipefail

EXP=$1
: "${L1_DIR:?}" "${LEGACY_DIR:?}"
OUT="${L1_DIR}/${EXP}"
IN="${L1_DIR}/inputs/${EXP}"

source "${LEGACY_DIR}/env.sh"
echo "[l1_task] ${EXP} task ${SLURM_ARRAY_TASK_ID} on $(hostname) at $(date -Is)"

case "$EXP" in
  endo_20250418)
    # dedup_map_folder.sh loops over a folder in one job; run its per-file step
    # (dedup_mapping.sh) for one sample, with the 15 threads the paper sbatch used.
    # Each sample gets its own base dir so the shared log files don't collide.
    FQ=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${IN}/fastqs.txt")
    S=$(basename "$FQ" .fq)
    bash "${LEGACY_DIR}/endogenous/20250418/dedup_mapping.sh" "$FQ" "$REF_ENDO" "${OUT}/${S}" 15
    ;;
  endo_20251022)
    bash "${LEGACY_DIR}/endogenous/20251022/trim_map_dedup.sh" \
      --map_file "${IN}/sample_map.txt" --input_dir "$FQ_1022" \
      --output_dir "$OUT" --ref_fasta "$REF_ENDO"
    ;;
  mpra_invitro)
    # input_dir holds the paper's concatenated FASTQs (find-order concatenation, see
    # dev/legacy/SOURCES.md); no barcode subdir exists there, so the script uses them as is.
    bash "${LEGACY_DIR}/mpra_invitro/trim_map_dedup_mpra.sh" \
      --map_file "${IN}/sample_map.txt" --input_dir "$FQ_INVITRO" \
      --output_dir "$OUT" --ref_fasta "$REF_POOL1"
    ;;
  mpra_incell)
    # The script indexes the sorted *.fq files of input_dir by task ID; input_dir is a
    # staging dir of symlinks to the selected samples only.
    bash "${LEGACY_DIR}/mpra_incell/trim_map_mpra.sh" \
      --input_dir "${IN}/fq" --output_dir "$OUT" --ref_fasta "$REF_POOL1"
    ;;
  *) echo "unknown experiment $EXP" >&2; exit 1 ;;
esac
echo "[l1_task] done $(date -Is)"
