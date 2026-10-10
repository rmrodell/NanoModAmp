#!/bin/bash
# Golden plan §4.6 step 1, check V1: run the legacy chain (preprocess → count → sample_name)
# on the FULL selection (select/<exp>/full, before the G1-d cap) with the full reference and
# the selected-chr BED, then compare with the paper counts (v1_compare.py, rule in
# docs/DECISIONS.md G1-c 2026-10-06).
#
# Usage: sbatch -p normal -c 8 --mem=32G -t 04:00:00 run_v1.sh
# Writes only under $GOLDEN_DIR/v1/.
set -euo pipefail
HERE=${HERE:-/home/users/rodell/NanoModAmp/dev/golden}
LEGACY_DIR=$(cd "${HERE}/../legacy" && pwd)
GOLDEN_DIR=${GOLDEN_DIR:-/scratch/users/rodell/NanoModAmp_build/golden}
SEL=${GOLDEN_DIR}/select
L1IN=${GOLDEN_DIR}/l1/inputs
V1=${GOLDEN_DIR}/v1
THREADS=${SLURM_CPUS_PER_TASK:-1}
REF_ENDO=/home/groups/nicolemm/rodell/fastas/endoG1G2.fasta
REF_POOL1=/home/groups/nicolemm/rodell/pool1/pool1_cleaned_noadapters.fasta
source ${HERE}/stage_refs.sh   # fresh .fai copies (pool1 .fai is stale)
BED_ENDO=/home/groups/nicolemm/rodell/fastas/endoG1G2_delpos.bed
BED_POOL1=/home/groups/nicolemm/rodell/pool1/pool1_cleaned_delpos_noadapters.bed
CHR_MPRA=/scratch/users/rodell/NanoModAmp_build/pool1_chr.txt
mkdir -p ${V1}/{inputs,beds}

# Selected-chr BEDs, coordinates unchanged (§4.7)
awk 'NR==FNR{k[$1]=1; next} ($1 in k)' ${L1IN}/endo_20250418/chrs.txt $BED_ENDO > ${V1}/beds/endo_20250418.bed
awk 'NR==FNR{k[$1]=1; next} ($1 in k)' ${L1IN}/endo_20251022/chrs.txt $BED_ENDO > ${V1}/beds/endo_20251022.bed
awk 'NR==FNR{if($1!="")k[$1]=1; next} ($1 in k)' $CHR_MPRA $BED_POOL1 > ${V1}/beds/mpra.bed

# Run a legacy array-task script once per sample-map line (it only uses the task ID as a line index).
run_tasks() {  # n_tasks script args...
  local n=$1; shift
  for i in $(seq 1 "$n"); do SLURM_ARRAY_TASK_ID=$i SLURM_CPUS_PER_TASK=$THREADS bash "$@"; done
}

(
source ${LEGACY_DIR}/env.sh

# --- endogenous 20250418: per-file script expects <sample>.fq; counts on pre-dedup sorted/ -----
E=endo_20250418; mkdir -p ${V1}/inputs/$E ${V1}/$E/count_bams
for gz in ${SEL}/$E/full/*.fastq.gz; do s=$(basename "$gz" .fastq.gz)
  zcat "$gz" > ${V1}/inputs/$E/$s.fq
  bash ${LEGACY_DIR}/endogenous/20250418/dedup_mapping.sh ${V1}/inputs/$E/$s.fq $REF_ENDO ${V1}/$E/$s $THREADS
  ln -sf ${V1}/$E/$s/minimap2/sorted/${s}_UMI_sort.bam ${V1}/$E/$s/minimap2/sorted/${s}_UMI_sort.bam.bai ${V1}/$E/count_bams/
done
mkdir -p ${V1}/$E/counts
bash ${LEGACY_DIR}/endogenous/20250418/DelDetect.sh ${V1}/$E/count_bams ${V1}/$E/counts $REF_ENDO ${V1}/beds/$E.bed

# --- endogenous 20251022: trim_map_dedup.sh finds <sample>.fastq.gz in input_dir --------------
E=endo_20251022
run_tasks "$(wc -l < ${L1IN}/$E/sample_map.txt)" ${LEGACY_DIR}/endogenous/20251022/trim_map_dedup.sh \
  --map_file ${L1IN}/$E/sample_map.txt --input_dir ${SEL}/$E/full --output_dir ${V1}/$E --ref_fasta $REF_ENDO
bash ${LEGACY_DIR}/endogenous/20251022/BIDdetect.sh -b ${V1}/$E/deduplicated_bam -o ${V1}/$E/counts \
  -r $REF_ENDO -e ${V1}/beds/$E.bed -n celltype_vector_rep_treat

# --- MPRA in vitro: no barcode subdirs in input_dir, so <sample>.fastq.gz is used directly -----
E=mpra_invitro
run_tasks "$(wc -l < ${L1IN}/$E/sample_map.txt)" ${LEGACY_DIR}/mpra_invitro/trim_map_dedup_mpra.sh \
  --map_file ${L1IN}/$E/sample_map.txt --input_dir ${SEL}/$E/full --output_dir ${V1}/$E --ref_fasta $REF_POOL1
bash ${LEGACY_DIR}/mpra_invitro/BIDdetect.sh -b ${V1}/$E/deduplicated_bam -o ${V1}/$E/counts \
  -r $REF_POOL1 -e ${V1}/beds/mpra.bed -n celltype_vector_treat_rep

# --- MPRA in cellulo: script indexes sorted *.fq of input_dir ---------------------------------
E=mpra_incell; mkdir -p ${V1}/inputs/$E
for gz in ${SEL}/$E/full/*.fastq.gz; do zcat "$gz" > ${V1}/inputs/$E/$(basename "$gz" .fastq.gz).fq; done
run_tasks "$(ls ${V1}/inputs/$E/*.fq | wc -l)" ${LEGACY_DIR}/mpra_incell/trim_map_mpra.sh \
  --input_dir ${V1}/inputs/$E --output_dir ${V1}/$E --ref_fasta $REF_POOL1
bash ${LEGACY_DIR}/mpra_incell/BIDdetect.sh -b ${V1}/$E/final_bam -o ${V1}/$E/counts \
  -r $REF_POOL1 -e ${V1}/beds/mpra.bed -n celltype_vector_rep_treat
) > ${V1}/legacy_run.log 2>&1

# Decompressed .fq copies were only needed as legacy-script inputs
rm -f ${V1}/inputs/endo_20250418/*.fq ${V1}/inputs/mpra_incell/*.fq

python3 ${HERE}/v1_compare.py --v1 ${V1} --out ${V1}/V1_results.tsv
