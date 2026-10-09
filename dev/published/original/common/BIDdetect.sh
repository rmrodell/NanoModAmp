#!/bin/bash

#
# BIDdetect Pipeline
#
# Description:
#   This script automates the process of running the BIDdetect pipeline. It iterates
#   through sorted BAM files in a source directory, runs 'bam_counts_fast.R' for each,
#   aggregates the results, and then runs 'sample_name.R' for final processing.
#

# Stop the script if any command fails or if an unset variable is used
set -eu -o pipefail

# Get the directory where this script is located to find the R scripts
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# --- 1. Argument Parsing and Usage ---

# Function to display usage information
usage() {
    echo "Usage: $0 -b <bam_dir> -o <output_dir> -r <ref_fasta> -e <bed_file> [OPTIONS]"
    echo ""
    echo "Required Arguments:"
    echo "  -b, --bam_dir <dir>      Directory containing sorted and indexed BAM files."
    echo "  -o, --output_dir <dir>   Directory where all outputs will be created."
    echo "  -r, --ref_fasta <file>   Path to the reference FASTA file."
    echo "  -e, --bed_file <file>    Path to the BED file with regions of interest."
    echo ""
    echo "Optional Arguments:"
    echo "  -n, --col_names <str>    Column names string for the final R script."
    echo "                           (Default: 'celltype_vector_rep_treat')"
    echo "  -h, --help               Display this help message and exit."
    echo ""
    echo "Example:"
    echo "  $0 --bam_dir ./bams/ --output_dir ./results/ --ref_fasta ref.fa --bed_file regions.bed"
}

# --- Set Default Optional Arguments ---
COLUMN_NAMES="celltype_vector_rep_treat"

# --- Parse Arguments using getopt ---
SHORT_OPTS="b:o:r:e:n:h"
LONG_OPTS="bam_dir:,output_dir:,ref_fasta:,bed_file:,col_names:,help"
PARSED_OPTS=$(getopt -o "$SHORT_OPTS" --long "$LONG_OPTS" -n "$0" -- "$@")
if [ $? -ne 0 ]; then usage; exit 1; fi
eval set -- "$PARSED_OPTS"

while true; do
    case "$1" in
        -b|--bam_dir)    BAM_DIR="$2"; shift 2 ;;
        -o|--output_dir) DEST_DIR="$2"; shift 2 ;;
        -r|--ref_fasta)  REF_FA="$2"; shift 2 ;;
        -e|--bed_file)   BED_FILE="$2"; shift 2 ;;
        -n|--col_names)  COLUMN_NAMES="$2"; shift 2 ;;
        -h|--help)       usage; exit 0 ;;
        --)              shift; break ;;
        *)               echo "Internal error!"; exit 1 ;;
    esac
done

# --- 2. Input Validation ---

# Check for missing required arguments and validate paths
if [ -z "${BAM_DIR-}" ] || [ -z "${DEST_DIR-}" ] || [ -z "${REF_FA-}" ] || [ -z "${BED_FILE-}" ]; then
    echo "Error: One or more required arguments are missing."
    usage
    exit 1
fi
if [ ! -d "$BAM_DIR" ]; then echo "Error: BAM directory not found: $BAM_DIR"; exit 1; fi
if [ ! -f "$REF_FA" ]; then echo "Error: Reference FASTA not found: $REF_FA"; exit 1; fi
if [ ! -f "$BED_FILE" ]; then echo "Error: BED file not found: $BED_FILE"; exit 1; fi

# --- 3. Setup Environment and Logging ---

# Load required modules
ml R

# Create output directories
LOG_DIR="${DEST_DIR}/logs"
INTERMEDIATE_DIR="${DEST_DIR}/intermediate_counts"
mkdir -p "$DEST_DIR" "$LOG_DIR" "$INTERMEDIATE_DIR"

# Set up comprehensive logging to a timestamped file and the console
LOG_FILE="${LOG_DIR}/BIDdetect_$(date +%Y%m%d_%H%M%S).log"
exec &> >(tee -a "$LOG_FILE")

echo "--- BIDdetect Pipeline Started: $(date) ---"
echo "BAM Directory:        $BAM_DIR"
echo "Output Directory:     $DEST_DIR"
echo "Reference FASTA:      $REF_FA"
echo "BED File:             $BED_FILE"
echo "Final Column Names:   $COLUMN_NAMES"
echo "Log file:             $LOG_FILE"
echo "-------------------------------------------------"

# --- 4. Main Processing Loop ---

# Define the master file for aggregated counts
master_count_file="${DEST_DIR}/BIDdetect_counts.txt"

# Initialize the master file with headers if it doesn't exist
if [ ! -f "$master_count_file" ]; then
    echo "Creating new master count file: $master_count_file"
    echo -e "sample\tchr\tpos\tgene\ttotalReads\tA.count\tC.count\tG.count\tT.count\tDeletion.count\tInsertion.count\tref\tkmer\tstrand\tdelrate" > "$master_count_file"
fi

# Check if there are any BAM files to process
if ! ls "${BAM_DIR}"/*.bam 1> /dev/null 2>&1; then
    echo "Error: No .bam files found in '$BAM_DIR'."
    exit 1
fi

# Loop through each .bam file in the input directory
for bam_file in "$BAM_DIR"/*.bam; do

    # Derives sample name from file, removing '_sort.bam' or just '.bam'
    sample=$(basename "${bam_file}" .bam | sed 's/_sort$//')
    echo "--- Processing sample: ${sample} ---"

    # Define path for the intermediate result file for this sample
    intermediate_result_file="${INTERMEDIATE_DIR}/${sample}_counts.txt"

    echo "Running R script to generate counts..."
    Rscript "${SCRIPT_DIR}/bam_counts_fast.R" \
        --bedFile "$BED_FILE" \
        --bamFile "$bam_file" \
        --referenceFasta "$REF_FA" \
        --outputFile "$intermediate_result_file"

    # Check if the intermediate file was created and has content before appending
    if [ -s "$intermediate_result_file" ]; then
        echo "Appending results to master count file..."
        # Append the results, skipping the header, and add the sample name
        awk -v sample="$sample" 'NR > 1 {print sample "\t" $0}' "$intermediate_result_file" >> "$master_count_file"
        echo "Finished processing sample: ${sample}"
    else
        # If the file is empty or missing, log a warning and continue
        echo "WARNING: No counts generated for sample ${sample}. Skipping append step."
    fi
    
done

echo "-------------------------------------------------"
echo "All samples have been processed."

# --- 5. Final Processing and Cleanup ---

final_output_file="${DEST_DIR}/BIDdetect_data.txt"

echo "Running final R script to format sample names..."
Rscript "${SCRIPT_DIR}/sample_name.R" "$master_count_file" "$final_output_file" "$COLUMN_NAMES"

# Optional: Uncomment the line below to automatically clean up intermediate files
# rm -rf "$INTERMEDIATE_DIR"

echo "--- Pipeline finished successfully! ---"
echo "Final formatted data is at: $final_output_file"
echo "Aggregated raw counts are at: $master_count_file"



