#!/bin/bash
# Golden package: rerun of the PUBLISHED pipeline (PUS7regulation2026 @ 4fd285b, dev/published/)
# on the final (capped) package FASTQs, saved next to the legacy rerun for comparison with the
# new pipeline (Becca, 2026-10-06). Not binding: the binding expected output is the legacy
# rerun (V2).
#
#   endogenous    Figure2 trim_map_dedup.sh (20251022) + CONSTRUCTED 3'-only variant (20250418, D31)
#                 (both runs, one experiment) → BIDdetect → analysis_endo.R
#   mpra_invitro  Figure3 trim_map_dedup_mpra.sh → BIDdetect → modification_analysis.R
#   mpra_incell   CONSTRUCTED trim_map_mpra_published.sh (pool adapters, no UMI/dedup; see
#                 dev/published/SOURCES.md) → BIDdetect → incell_analysis.R (cell-type analysis off)
#
# Endogenous inputs are staged under normalised names so the two runs form one table with the
# vector labels analysis_endo.R expects (WT/KD) and run-unique reps:
#   20250418: P36, pLKO → WT; P97, shPUS7 → KD; rep r → 0418r<r>.   20251022: rep r → 1022r<r>.
#
# Usage: sbatch -p normal -c 8 --mem=32G -t 04:00:00 run_published.sh
# Env:   SELECT_DIR (capped FASTQs, default $GOLDEN_DIR/select), OUT_DIR, WORK_DIR, EXPERIMENTS
# Writes $OUT_DIR/<exp>/{read_funnel.tsv,BIDdetect_data.txt,site_calls/,plots_png/,bam_md5.txt}.
# Idempotent: each experiment's work and output dirs are rebuilt from scratch.
set -euo pipefail
HERE=${HERE:-/home/users/rodell/NanoModAmp/dev/golden}
PUB=$(cd "${HERE}/../published" && pwd)
GOLDEN_DIR=${GOLDEN_DIR:-/scratch/users/rodell/NanoModAmp_build/golden}
SELECT_DIR=${SELECT_DIR:-${GOLDEN_DIR}/select}
OUT_DIR=${OUT_DIR:-${GOLDEN_DIR}/expected/published_rerun}
WORK_DIR=${WORK_DIR:-${GOLDEN_DIR}/published_work}
EXPERIMENTS=${EXPERIMENTS:-"endogenous mpra_invitro mpra_incell"}
THREADS=${SLURM_CPUS_PER_TASK:-1}
L1IN=${GOLDEN_DIR}/l1/inputs
REF_ENDO=/home/groups/nicolemm/rodell/fastas/endoG1G2.fasta
REF_POOL1=/home/groups/nicolemm/rodell/pool1/pool1_cleaned_noadapters.fasta
BED_ENDO=/home/groups/nicolemm/rodell/fastas/endoG1G2_delpos.bed
BED_POOL1=/home/groups/nicolemm/rodell/pool1/pool1_cleaned_delpos_noadapters.bed
CHR_ENDO=/scratch/users/rodell/NanoModAmp_build/endoBID_data/endo_chr.txt
CHR_MPRA=/scratch/users/rodell/NanoModAmp_build/pool1_chr.txt
RSCRIPT43=/share/software/user/open/R/4.3.2/bin/Rscript

# Run a published array-task script once per sample-map line (it only uses the task ID as a line index).
run_tasks() {  # map script args...
  local map=$1; shift
  local n i; n=$(wc -l < "$map")
  for i in $(seq 1 "$n"); do SLURM_ARRAY_TASK_ID=$i SLURM_CPUS_PER_TASK=$THREADS bash "$@" --map_file "$map"; done
}

# Selected-chr BED, coordinates unchanged
sel_bed() { awk 'NR==FNR{if($1!="")k[$1]=1; next} ($1 in k)' <(cat "$1"; echo) "$2"; }

# read funnel (all samples) from the per-sample run_metrics tables
funnel() {  # prep_dir out
  local f first=1
  for f in "$1"/*/reports/*.run_metrics.tsv; do
    if [ $first = 1 ]; then head -1 "$f"; first=0; fi
    tail -n +2 "$f"
  done > "$2"
}

collect() {  # exp counts_dir bam_dir site_calls_src
  local exp=$1 out=${OUT_DIR}/$1
  cp "$2/BIDdetect_data.txt" "$out/"
  (cd "$3" && md5sum *.bam) > "$out/bam_md5.txt"
  mkdir -p "$out/site_calls" "$out/plots_png"
  (cd "$4" && find . -type f -name '*.png' | sort | while read -r f; do cp "$f" "$out/plots_png/$(echo "${f#./}" | tr / _)"; done)
  (cd "$4" && find . -type f \( -name '*.tsv' -o -name '*.csv' -o -name '*.txt' \) | sort | while read -r f; do
     cp "$f" "$out/site_calls/$(echo "${f#./}" | tr / _)"; done)
}

# References are staged with a freshly built .fai: the shipped pool1 .fai (2025-11-07) predates a
# rewrite of the FASTA (2026-01-26), so its offsets are wrong and Rsamtools FaFile reads garbage.
stage_ref() {  # src dest_dir -> prints staged path
  mkdir -p "$2"; cp "$1" "$2/"
  samtools faidx "$2/$(basename "$1")"
  echo "$2/$(basename "$1")"
}
set +u; module load biology samtools/1.16.1 >/dev/null 2>&1; set -u

run_exp() {
  local exp=$1 W O REF
  W=${WORK_DIR}/$exp; O=${OUT_DIR}/$exp
  rm -rf "$W" "$O"; mkdir -p "$W"/{fastq,beds} "$O"
  echo "=== ${exp} $(date -Is)"
  case $exp in endogenous) REF=$(stage_ref "$REF_ENDO" "$W/ref") ;; *) REF=$(stage_ref "$REF_POOL1" "$W/ref") ;; esac
  case $exp in
    endogenous)
      : > "$W/sample_map_20250418.txt"; : > "$W/sample_map_20251022.txt"
      for run in 20250418 20251022; do
        for gz in ${SELECT_DIR}/endo_${run}/capped/*.fastq.gz; do
          IFS=_ read -r ct vec rep trt <<< "$(basename "$gz" .fastq.gz)"
          case $vec in P36|pLKO) vec=WT ;; P97|shPUS7) vec=KD ;; esac
          new="${ct}_${vec}_${run:4}r${rep}_${trt}"
          ln -s "$gz" "$W/fastq/${new}.fastq.gz"
          echo "${new}:none" >> "$W/sample_map_${run}.txt"
          echo -e "${run}\t$(basename "$gz" .fastq.gz)\t${new}" >> "$W/name_map.tsv"
        done
      done
      sel_bed "$CHR_ENDO" "$BED_ENDO" > "$W/beds/targets.bed"
      ( set +u; source ${PUB}/env.sh
        # 20250418: constructed 3'-only variant (D31/R-29; reads lack the 5' ONT adapter)
        run_tasks "$W/sample_map_20250418.txt" ${PUB}/endogenous/trim_map_dedup_3prime_only.sh \
          --input_dir "$W/fastq" --output_dir "$W/prep" --ref_fasta "$REF"
        run_tasks "$W/sample_map_20251022.txt" ${PUB}/endogenous/trim_map_dedup.sh \
          --input_dir "$W/fastq" --output_dir "$W/prep" --ref_fasta "$REF"
        bash ${PUB}/common/BIDdetect.sh -b "$W/prep/deduplicated_bam" -o "$W/counts" \
          -r "$REF" -e "$W/beds/targets.bed" -n celltype_vector_rep_treat )
      mkdir -p "$W/site_calls"
      ( set +u; unset R_HOME; module load R/4.3.2 >/dev/null 2>&1  # swap, not purge: purge drops libgfortran (gcc/12)
        cd "$W/site_calls" && PUB_ENDO_DIR="$W/counts" $RSCRIPT43 ${PUB}/endogenous/analysis_endo.R )
      cp "$W/name_map.tsv" "$O/"
      funnel "$W/prep" "$O/read_funnel.tsv"
      collect $exp "$W/counts" "$W/prep/deduplicated_bam" "$W/site_calls"
      ;;
    mpra_invitro)
      for gz in ${SELECT_DIR}/$exp/capped/*.fastq.gz; do ln -s "$gz" "$W/fastq/"; done
      sel_bed "$CHR_MPRA" "$BED_POOL1" > "$W/beds/targets.bed"
      ( set +u; source ${PUB}/env.sh
        run_tasks "${L1IN}/mpra_invitro/sample_map.txt" ${PUB}/mpra_invitro/trim_map_dedup_mpra.sh \
          --input_dir "$W/fastq" --output_dir "$W/prep" --ref_fasta "$REF"
        bash ${PUB}/common/BIDdetect.sh -b "$W/prep/deduplicated_bam" -o "$W/counts" \
          -r "$REF" -e "$W/beds/targets.bed" -n celltype_vector_treat_rep )
      ( set +u; unset R_HOME; module load R/4.3.2 >/dev/null 2>&1  # swap, not purge: purge drops libgfortran (gcc/12)
        $RSCRIPT43 ${PUB}/mpra_invitro/modification_analysis.R --input "$W/counts/BIDdetect_data.txt" \
          --outdir "$W/site_calls" --prefix invitro_delpos --cores "$THREADS" --plot_all_sites )
      funnel "$W/prep" "$O/read_funnel.tsv"
      collect $exp "$W/counts" "$W/prep/deduplicated_bam" "$W/site_calls"
      ;;
    mpra_incell)
      : > "$W/sample_map.txt"
      for gz in ${SELECT_DIR}/$exp/capped/*.fastq.gz; do
        ln -s "$gz" "$W/fastq/"; echo "$(basename "$gz" .fastq.gz):none" >> "$W/sample_map.txt"
      done
      sel_bed "$CHR_MPRA" "$BED_POOL1" > "$W/beds/targets.bed"
      ( set +u; source ${PUB}/env.sh
        run_tasks "$W/sample_map.txt" ${PUB}/mpra_incell/trim_map_mpra_published.sh \
          --input_dir "$W/fastq" --output_dir "$W/prep" --ref_fasta "$REF"
        bash ${PUB}/common/BIDdetect.sh -b "$W/prep/final_bam" -o "$W/counts" \
          -r "$REF" -e "$W/beds/targets.bed" -n celltype_vector_rep_treat )
      mkdir -p "$W/site_calls"
      ( set +u; unset R_HOME; module load R/4.3.2 >/dev/null 2>&1  # swap, not purge: purge drops libgfortran (gcc/12)
        cd "$W/site_calls" && PUB_INPUT="$W/counts/BIDdetect_data.txt" PUB_OUTDIR="$W/site_calls/incell" \
          PUB_CORES="$THREADS" PUB_RUN_CELLTYPE=FALSE $RSCRIPT43 ${PUB}/mpra_incell/incell_analysis.R )
      funnel "$W/prep" "$O/read_funnel.tsv"
      collect $exp "$W/counts" "$W/prep/final_bam" "$W/site_calls"
      ;;
  esac
}

# Each experiment runs in its own subshell with set -e, so one failure does not stop the others;
# the outcome is recorded in $OUT_DIR/status.tsv and the job exits non-zero if any failed.
mkdir -p "$OUT_DIR"; : > "$OUT_DIR/status.tsv"; fail=0
for exp in $EXPERIMENTS; do
  set +e; ( set -e; run_exp "$exp" ); rc=$?; set -e
  echo -e "${exp}\t$([ $rc = 0 ] && echo ok || echo "FAILED (exit $rc)")" | tee -a "$OUT_DIR/status.tsv"
  [ $rc = 0 ] || fail=1
done
echo "=== done $(date -Is)"
exit $fail
