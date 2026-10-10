#!/bin/bash
# Golden plan §5 check Repro (selection/extraction stage): rebuild every sample's FASTQs from the
# same inputs and seed into select_repro/ and require byte-identical files to select/.
# The package-level check (MANIFEST.md5) is repeated by assemble_package.py.
#
# Usage: sbatch -p normal -c 8 --mem=32G -t 02:00:00 run_repro.sh   (after run_select.sh)
set -euo pipefail
HERE=${HERE:-/home/users/rodell/NanoModAmp/dev/golden}
GOLDEN_DIR=${GOLDEN_DIR:-/scratch/users/rodell/NanoModAmp_build/golden}
SEED=20261004 CAP=250 BG=150
TASKS=${GOLDEN_DIR}/select/tasks.tsv
OUT=${GOLDEN_DIR}/select_repro
module load biology samtools/1.16.1 >/dev/null 2>&1
rm -rf "$OUT"; mkdir -p "$OUT"

export HERE OUT SEED CAP BG
run_one() {
  IFS=$'\t' read -r exp s bam fq chrs umi <<< "$1"
  python3 "${HERE}/select_extract.py" --sample "$s" --bam "$bam" --fastq "$fq" --chrs "$chrs" \
    --outdir "${OUT}/${exp}" --seed "$SEED" --cap "$CAP" --background "$BG" --umi-suffix "$umi" > /dev/null
}
export -f run_one
tr '\n' '\0' < "$TASKS" | xargs -0 -P "${SLURM_CPUS_PER_TASK:-1}" -I{} bash -c 'run_one "$1"' _ {}

( cd "${GOLDEN_DIR}/select" && find . -path '*/full/*' -o -path '*/capped/*' -o -path '*/ids/*.txt' | sort | xargs md5sum ) > "${OUT}/orig.md5"
( cd "$OUT" && find . -path '*/full/*' -o -path '*/capped/*' -o -path '*/ids/*.txt' | sort | xargs md5sum ) > "${OUT}/repro.md5"
if cmp -s "${OUT}/orig.md5" "${OUT}/repro.md5"; then
  echo "Repro PASS: $(wc -l < "${OUT}/repro.md5") files byte-identical"
else
  echo "Repro FAIL:"; diff "${OUT}/orig.md5" "${OUT}/repro.md5" | head -20; exit 1
fi
