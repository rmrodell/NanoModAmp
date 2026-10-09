#!/bin/bash
set -e
set -u

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Get the root directory of the git repo (assuming the script is in the root)
REPO_ROOT="$SCRIPT_DIR"

# This trims one set of adapters

# Check if required arguments are provided
if [ $# -lt 3 ]; then
    echo "Usage: $0 <input_directory> <reference_fasta> <base_directory> [threads]"
    echo "Example: $0 /path/to/fastq/dir /path/to/reference.fasta /scratch/users/username/output_dir 7"
    echo "Note: If threads not specified, default value of 1 will be used"
    exit 1
fi

# Assign command line arguments to variables
INPUT_DIR=$1
REFERENCE_FASTA=$2
BASE_DIR=$3
THREADS=${4:-1}

# Verify input directory exists
if [ ! -d "${INPUT_DIR}" ]; then
    echo "Error: Input directory ${INPUT_DIR} not found"
    exit 1
fi

# Verify reference exists
if [ ! -f "${REFERENCE_FASTA}" ]; then
    echo "Error: Reference fasta ${REFERENCE_FASTA} not found"
    exit 1
fi

# Process each .fq file in the input directory
for fastq in ${INPUT_DIR}/*.fq; do
    if [ -f "$fastq" ]; then
        echo "Processing $fastq"
        bash "${REPO_ROOT}/dedup_mapping.sh" "$fastq" "$REFERENCE_FASTA" "$BASE_DIR" "$THREADS"
    fi
done

# Create a combined text summary file
echo "Creating combined text summary..."
TEXT_SUMMARY="${BASE_DIR}/logs/combined_summary.txt"
echo "Combined Processing Summary" > ${TEXT_SUMMARY}
echo "----------------------------------------" >> ${TEXT_SUMMARY}

# Combine all individual summaries
for summary in ${BASE_DIR}/logs/*_summary.txt; do
    if [ -f "$summary" ]; then
        echo "" >> ${TEXT_SUMMARY}
        cat "$summary" >> ${TEXT_SUMMARY}
        echo "----------------------------------------" >> ${TEXT_SUMMARY}
    fi
done

# Create a combined summary table
echo "Creating combined summary table..."
TABLE_SUMMARY="${BASE_DIR}/logs/combined_summary_table.tsv"


echo
