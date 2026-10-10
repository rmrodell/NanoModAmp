#!/bin/bash
# Golden plan §4.3–4.5 for all approved samples: read-ID assignment from the L1 rerun
# pre-dedup BAMs, seeded background, extraction (full selection for V1, capped for the
# package). One job, samples run in parallel (they take 1–5 min each).
#
# Usage: sbatch -p normal -c 8 --mem=32G -t 02:00:00 run_select.sh
# Writes $GOLDEN_DIR/select/<experiment>/{full,capped,ids}/ and select_summary.tsv
set -euo pipefail
HERE=${HERE:-/home/users/rodell/NanoModAmp/dev/golden}
GOLDEN_DIR=${GOLDEN_DIR:-/scratch/users/rodell/NanoModAmp_build/golden}
L1=${GOLDEN_DIR}/l1
SEED=20261004 CAP=250 BG=150          # build settings (golden_inputs.yaml, G1-d 2026-10-06)
BID=/oak/stanford/groups/nicolemm/rodell/BIDamplicon
module load biology samtools/1.16.1 >/dev/null 2>&1

# experiment  sample  pre-dedup BAM  raw FASTQ (as read by the legacy script)  chrs  umi_len
TASKS=${GOLDEN_DIR}/select/tasks.tsv
mkdir -p ${GOLDEN_DIR}/select
: > "$TASKS"
while read -r fq; do s=$(basename "$fq" .fq)
  echo -e "endo_20250418\t$s\t${L1}/endo_20250418/$s/minimap2/sorted/${s}_UMI_sort.bam\t$fq\t${L1}/inputs/endo_20250418/chrs.txt\t10"
done < ${L1}/inputs/endo_20250418/fastqs.txt >> "$TASKS"
while IFS=: read -r s _; do
  echo -e "endo_20251022\t$s\t${L1}/endo_20251022/retained/${s}_mapped_sorted.bam\t${BID}/endo_CellTypeSpec/20251022_endoBID/fastq/$s.fastq.gz\t${L1}/inputs/endo_20251022/chrs.txt\t10"
done < ${L1}/inputs/endo_20251022/sample_map.txt >> "$TASKS"
while IFS=: read -r s _; do
  echo -e "mpra_invitro\t$s\t${L1}/mpra_invitro/retained/${s}_mapped_sorted.bam\t${BID}/Pool1/invitro/maydata/00_concatenated_fastqs/$s.fastq.gz\t/scratch/users/rodell/NanoModAmp_build/pool1_chr.txt\t10"
done < ${L1}/inputs/mpra_invitro/sample_map.txt >> "$TASKS"
for f in ${L1}/inputs/mpra_incell/fq/*.fq; do s=$(basename "$f" .fq)
  echo -e "mpra_incell\t$s\t${L1}/mpra_incell/final_bam/$s.bam\t$(readlink -f "$f")\t/scratch/users/rodell/NanoModAmp_build/pool1_chr.txt\t0"
done >> "$TASKS"
[ "$(wc -l < "$TASKS")" = 68 ] || { echo "expected 68 tasks" >&2; exit 1; }

export HERE GOLDEN_DIR SEED CAP BG
run_one() {
  IFS=$'\t' read -r exp s bam fq chrs umi <<< "$1"
  python3 "${HERE}/select_extract.py" --sample "$s" --bam "$bam" --fastq "$fq" --chrs "$chrs" \
    --outdir "${GOLDEN_DIR}/select/${exp}" --seed "$SEED" --cap "$CAP" --background "$BG" \
    --umi-suffix "$umi" > /dev/null
}
export -f run_one
tr '\n' '\0' < "$TASKS" | xargs -0 -P "${SLURM_CPUS_PER_TASK:-1}" -I{} bash -c 'run_one "$1"' _ {}

# Summary table (per sample and chr) and sizes
python3 - "$GOLDEN_DIR/select" <<'EOF'
import glob, json, os, sys
d = sys.argv[1]
with open(os.path.join(d, "select_summary.tsv"), "w") as out:
    out.write("experiment\tsample\tchr\tselected\tcapped\tbackground\traw_reads\tL2\n")
    tot = {}
    for f in sorted(glob.glob(os.path.join(d, "*", "ids", "*.summary.json"))):
        exp = f.split(os.sep)[-3]; j = json.load(open(f))
        ok = j["L2_full"]["wanted"] == j["L2_full"]["found"] and j["L2_capped"]["wanted"] == j["L2_capped"]["found"]
        for c, v in sorted(j["per_chr"].items()):
            out.write("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\n" % (exp, j["sample"], c, v["selected"], v["capped"], j["background_reads"], j["raw_reads"], "pass" if ok else "FAIL"))
for kind in ("full", "capped"):
    for exp in sorted(os.listdir(d)):
        fs = glob.glob(os.path.join(d, exp, kind, "*.fastq.gz"))
        if fs:
            print("%-8s %-14s %3d files %8.2f MB" % (kind, exp, len(fs), sum(os.path.getsize(f) for f in fs) / 1e6))
EOF
