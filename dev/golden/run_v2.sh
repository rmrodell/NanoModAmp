#!/bin/bash
# Golden plan §4.6 step 2 + §4.8 legacy_rerun (V2, binding expected output): size check, then
# the full legacy chain on the FINAL package (select/<run>/capped): preprocessing → counting →
# sample_name → legacy site calling. Also copies the V1 counts to expected/legacy_full_selected.
#
# Usage: sbatch -p normal -c 8 --mem=32G -t 06:00:00 run_v2.sh
# Writes only under $GOLDEN_DIR/{v2,expected}/ and $GOLDEN_DIR/size_report.txt.
set -euo pipefail
HERE=${HERE:-/home/users/rodell/NanoModAmp/dev/golden}
LEGACY_DIR=$(cd "${HERE}/../legacy" && pwd)
GOLDEN_DIR=${GOLDEN_DIR:-/scratch/users/rodell/NanoModAmp_build/golden}
SEL=${GOLDEN_DIR}/select
L1IN=${GOLDEN_DIR}/l1/inputs
V1=${GOLDEN_DIR}/v1
V2=${GOLDEN_DIR}/v2
EXP_DIR=${GOLDEN_DIR}/expected
THREADS=${SLURM_CPUS_PER_TASK:-1}
RSCRIPT43=/share/software/user/open/R/4.3.2/bin/Rscript   # site calling; counting stays on R/4.2.0
REF_ENDO=/home/groups/nicolemm/rodell/fastas/endoG1G2.fasta
REF_POOL1=/home/groups/nicolemm/rodell/pool1/pool1_cleaned_noadapters.fasta
source ${HERE}/stage_refs.sh   # fresh .fai copies (pool1 .fai is stale)
BEDS=${V1}/beds   # selected-chr delpos BEDs written by run_v1.sh (coordinates unchanged)

# --- §4.6 step 2: size check (stop before V2 if over budget) ---------------------------------
python3 ${HERE}/check_size.py --golden ${GOLDEN_DIR} --report ${GOLDEN_DIR}/size_report.txt

mkdir -p ${V2}/inputs
run_tasks() {  # n_tasks script args...  (legacy array scripts only use the task ID as a line index)
  local n=$1; shift
  for i in $(seq 1 "$n"); do SLURM_ARRAY_TASK_ID=$i SLURM_CPUS_PER_TASK=$THREADS bash "$@"; done
}

# --- legacy preprocessing + counting, exactly as run_v1.sh but on capped/ ---------------------
(
source ${LEGACY_DIR}/env.sh

E=endo_20250418; mkdir -p ${V2}/inputs/$E ${V2}/$E/count_bams
for gz in ${SEL}/$E/capped/*.fastq.gz; do s=$(basename "$gz" .fastq.gz)
  zcat "$gz" > ${V2}/inputs/$E/$s.fq
  bash ${LEGACY_DIR}/endogenous/20250418/dedup_mapping.sh ${V2}/inputs/$E/$s.fq $REF_ENDO ${V2}/$E/$s $THREADS
  ln -sf ${V2}/$E/$s/minimap2/sorted/${s}_UMI_sort.bam ${V2}/$E/$s/minimap2/sorted/${s}_UMI_sort.bam.bai ${V2}/$E/count_bams/
done
mkdir -p ${V2}/$E/counts
bash ${LEGACY_DIR}/endogenous/20250418/DelDetect.sh ${V2}/$E/count_bams ${V2}/$E/counts $REF_ENDO ${BEDS}/$E.bed

E=endo_20251022
run_tasks "$(wc -l < ${L1IN}/$E/sample_map.txt)" ${LEGACY_DIR}/endogenous/20251022/trim_map_dedup.sh \
  --map_file ${L1IN}/$E/sample_map.txt --input_dir ${SEL}/$E/capped --output_dir ${V2}/$E --ref_fasta $REF_ENDO
bash ${LEGACY_DIR}/endogenous/20251022/BIDdetect.sh -b ${V2}/$E/deduplicated_bam -o ${V2}/$E/counts \
  -r $REF_ENDO -e ${BEDS}/$E.bed -n celltype_vector_rep_treat

E=mpra_invitro
run_tasks "$(wc -l < ${L1IN}/$E/sample_map.txt)" ${LEGACY_DIR}/mpra_invitro/trim_map_dedup_mpra.sh \
  --map_file ${L1IN}/$E/sample_map.txt --input_dir ${SEL}/$E/capped --output_dir ${V2}/$E --ref_fasta $REF_POOL1
bash ${LEGACY_DIR}/mpra_invitro/BIDdetect.sh -b ${V2}/$E/deduplicated_bam -o ${V2}/$E/counts \
  -r $REF_POOL1 -e ${BEDS}/mpra.bed -n celltype_vector_treat_rep

E=mpra_incell; mkdir -p ${V2}/inputs/$E
for gz in ${SEL}/$E/capped/*.fastq.gz; do zcat "$gz" > ${V2}/inputs/$E/$(basename "$gz" .fastq.gz).fq; done
run_tasks "$(ls ${V2}/inputs/$E/*.fq | wc -l)" ${LEGACY_DIR}/mpra_incell/trim_map_mpra.sh \
  --input_dir ${V2}/inputs/$E --output_dir ${V2}/$E --ref_fasta $REF_POOL1
bash ${LEGACY_DIR}/mpra_incell/BIDdetect.sh -b ${V2}/$E/final_bam -o ${V2}/$E/counts \
  -r $REF_POOL1 -e ${BEDS}/mpra.bed -n celltype_vector_rep_treat
) > ${V2}/legacy_preprocess_count.log 2>&1
rm -f ${V2}/inputs/endo_20250418/*.fq ${V2}/inputs/mpra_incell/*.fq

# --- legacy site calling (R 4.3.2 + ~/R/x86_64-pc-linux-gnu-library/4.3) ----------------------
(
# Swap to R/4.3.2 (brings gcc/12 for libgfortran.so.5); a purge would drop it. R_HOME may be pinned to 4.2.
unset R_HOME; module load R/4.3.2 >/dev/null 2>&1
$RSCRIPT43 -e 'for (p in c("blme","car","dplyr","readr","ggplot2","ggrepel","lme4","DescTools","furrr","argparse","ggtext","forcats")) suppressPackageStartupMessages(library(p, character.only=TRUE))'
SC=${V2}/site_calls; mkdir -p ${SC}/{mpra_invitro,mpra_incell,endogenous}

# in vitro: noPUS → input is applied by the script itself (as in the paper)
$RSCRIPT43 ${LEGACY_DIR}/mpra_invitro/modification_analysis.R -i ${V2}/mpra_invitro/counts/BIDdetect_data.txt \
  -o ${SC}/mpra_invitro --prefix invitro --cores $THREADS --plot_all_sites

# in cellulo: WT_mod + PUS7_dep only (cell-type analyses out of scope, Becca 2026-10-06)
(cd ${SC}/mpra_incell && LEGACY_INPUT=${V2}/mpra_incell/counts/BIDdetect_data.txt LEGACY_OUTDIR=${SC}/mpra_incell \
  LEGACY_CORES=$THREADS LEGACY_RUN_CELLTYPE=FALSE $RSCRIPT43 ${LEGACY_DIR}/mpra_incell/incell_analysis.R)

# endogenous: Figure 2 analysis_endo.R reads one combined BIDdetect_data.txt from LEGACY_ENDO_DIR and
# writes to the working directory; the combined table is built by legacy_combine_endo.py.
python3 ${HERE}/legacy_combine_endo.py ${SC}/endogenous/BIDdetect_data.txt \
  20250418=${V2}/endo_20250418/counts/deldetect_factors.txt 20251022=${V2}/endo_20251022/counts/BIDdetect_data.txt
(cd ${SC}/endogenous && LEGACY_ENDO_DIR=${SC}/endogenous $RSCRIPT43 ${LEGACY_DIR}/endogenous/site_calling/analysis_endo.R)
) > ${V2}/legacy_site_calling.log 2>&1

# --- expected/legacy_rerun/<exp>/ and expected/legacy_full_selected/<exp>/ ---------------------
python3 - "${GOLDEN_DIR}" <<'EOF'
import csv, glob, hashlib, os, shutil, sys
g = sys.argv[1]; v1 = os.path.join(g, "v1"); v2 = os.path.join(g, "v2"); sc = os.path.join(v2, "site_calls")
out = os.path.join(g, "expected", "legacy_rerun"); full = os.path.join(g, "expected", "legacy_full_selected")
for d in (out, full):  # idempotent: rebuild from v1/ and v2/ every time
    shutil.rmtree(d, ignore_errors=True)
md5 = lambda f: hashlib.md5(open(f, "rb").read()).hexdigest()
# run -> (package experiment, count table, BAM globs (count source first), funnel source)
RUNS = {
  "endo_20250418": ("endogenous", "counts/deldetect_factors.txt", ["*/minimap2/sorted/*_UMI_sort.bam", "*/minimap2/dedup/*_UMI_dedup.bam"], "summary"),
  "endo_20251022": ("endogenous", "counts/BIDdetect_data.txt", ["deduplicated_bam/*.bam", "retained/*_mapped_sorted.bam"], "metrics"),
  "mpra_invitro": ("mpra_invitro", "counts/BIDdetect_data.txt", ["deduplicated_bam/*.bam", "retained/*_mapped_sorted.bam"], "metrics"),
  "mpra_incell": ("mpra_incell", "counts/BIDdetect_data.txt", ["final_bam/*.bam"], "metrics"),
}
for run, (exp, table, bams, funnel) in RUNS.items():
    d = os.path.join(out, exp); os.makedirs(d, exist_ok=True)
    tag = run if exp == "endogenous" else ""
    shutil.copy(os.path.join(v2, run, table), os.path.join(d, ("%s_" % tag if tag else "") + "BIDdetect_data.txt"))
    with open(os.path.join(d, "bam_md5.txt"), "a") as fh:
        for pat in bams:
            for b in sorted(glob.glob(os.path.join(v2, run, pat))):
                fh.write("%s  %s/%s\n" % (md5(b), run, os.path.relpath(b, os.path.join(v2, run))))
    rows = []
    if funnel == "metrics":
        for f in sorted(glob.glob(os.path.join(v2, run, "*", "reports", "*.run_metrics.tsv"))):
            rows += [r for r in csv.reader(open(f), delimiter="\t")][1:]
    else:  # dedup_mapping.sh writes logs/<sample>_summary.txt ("Label: N" lines)
        for f in sorted(glob.glob(os.path.join(v2, run, "*", "logs", "*_summary.txt"))):
            s = os.path.basename(f)[:-len("_summary.txt")]
            for line in open(f):  # count lines only; percentage lines end in "%" and are skipped
                k, _, val = line.partition(":")
                if val.strip().isdigit():
                    rows.append([k.strip(), s, "", val.strip(), ""])
    with open(os.path.join(d, "read_funnel.tsv"), "a") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        if fh.tell() == 0:
            w.writerow(["run", "step", "sample", "file", "records", "size"])
        for r in rows:
            w.writerow([run] + r)
    # V1 counts (non-downsampled selection), counts only
    fd = os.path.join(full, exp); os.makedirs(fd, exist_ok=True)
    shutil.copy(os.path.join(v1, run, table), os.path.join(fd, ("%s_" % tag if tag else "") + "BIDdetect_data.txt"))
# site calls: tables + PNG only; list PDFs that have no PNG twin
for exp in ("mpra_invitro", "mpra_incell", "endogenous"):
    src = os.path.join(sc, exp); d = os.path.join(out, exp)
    os.makedirs(os.path.join(d, "site_calls"), exist_ok=True); os.makedirs(os.path.join(d, "plots_png"), exist_ok=True)
    pdf_only = []
    for root, _, files in os.walk(src):
        for f in files:
            p = os.path.join(root, f)
            if f.endswith((".tsv", ".csv", ".txt")) and f != "BIDdetect_data.txt":
                shutil.copy(p, os.path.join(d, "site_calls", f))
            elif f.endswith(".png"):
                shutil.copy(p, os.path.join(d, "plots_png", f))
            elif f.endswith(".pdf") and f != "Rplots.pdf" and not os.path.exists(p[:-4] + ".png"):
                pdf_only.append(os.path.relpath(p, src))
    open(os.path.join(d, "plots_pdf_only.txt"), "w").write("".join(x + "\n" for x in sorted(pdf_only)))
    if exp == "endogenous":
        shutil.copy(os.path.join(src, "BIDdetect_data.txt"), os.path.join(d, "combined_BIDdetect_data.txt"))
print("expected outputs written under", out)
EOF
du -sh ${EXP_DIR}/legacy_rerun/* ${EXP_DIR}/legacy_full_selected/*
