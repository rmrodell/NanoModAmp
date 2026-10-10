#!/bin/bash
# Golden plan §4.7 + §4.10: assemble tests/data/golden/ after the V2 legacy rerun and the
# published rerun have finished. Strict build into the repo; if inputs are missing it fails
# there and writes an --allow-missing build to $GOLDEN_DIR/assemble_partial/ for inspection.
#
# Usage: sbatch -p normal -c 1 --mem=4G -t 00:30:00 run_assemble.sh
set -uo pipefail
HERE=${HERE:-/home/users/rodell/NanoModAmp/dev/golden}
GOLDEN_DIR=${GOLDEN_DIR:-/scratch/users/rodell/NanoModAmp_build/golden}
module load biology samtools/1.16.1 htslib >/dev/null 2>&1

if python3 "${HERE}/assemble_package.py" --golden "$GOLDEN_DIR"; then
  echo "Package assembled in $(cd "${HERE}/../.." && pwd)/tests/data/golden"
else
  rc=$?
  echo "Strict assembly failed (exit $rc); writing partial build to ${GOLDEN_DIR}/assemble_partial" >&2
  python3 "${HERE}/assemble_package.py" --golden "$GOLDEN_DIR" --out "${GOLDEN_DIR}/assemble_partial" --allow-missing
  exit $rc
fi
