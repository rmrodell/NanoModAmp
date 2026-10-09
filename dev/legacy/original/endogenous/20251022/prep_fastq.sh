#!/bin/bash

# ---
# Concatenate Demultiplexed FASTQ Files
#
# Description:
#   This script processes demultiplexed nanopore data. For each sample
#   specified in a mapping file, it finds the corresponding source directory,
#   concatenates all gzipped FASTQ files within it, and saves the
#   resulting uncompressed FASTQ file to a destination directory.
#
# Dependencies:
#   - bash
#   - coreutils (cat, mkdir, tee, echo)
#   - gzip
#
# Input File Format (colon-separated, no header):
#   <output_filename_prefix>:<barcode_id>
#   HepG2_WT_1_BS:barcode75
#   293T_KD_1_input:barcode76
#
# Usage:
#   ./prep_fastq.sh -s /path/to/source_dir -d /path/to/dest_dir -b mapping_file.txt
#
# ---

# Function to display usage information
usage() {
    echo "Usage: $0 -s <source_dir> -d <destination_dir> -b <barcode_file>"
    echo "  -s  Source directory containing barcode subdirectories."
    echo "  -d  Destination directory for concatenated .fq files."
    echo "  -b  File mapping output filenames to barcodes (format: <filename>:<barcode>)."
    echo "  -h  Display this help message."
    exit 1
}

# Parse command-line options
while getopts ":s:d:b:h" opt; do
    case ${opt} in
        s) SOURCE_DIR="$OPTARG" ;;
        d) DEST_DIR="$OPTARG" ;;
        b) BARCODE_FILE="$OPTARG" ;;
        h) usage ;;
        \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
        :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done

# Check if all required arguments were provided
if [ -z "$SOURCE_DIR" ] || [ -z "$DEST_DIR" ] || [ -z "$BARCODE_FILE" ]; then
    echo "Error: Missing required arguments."
    usage
fi

# Check if the barcode mapping file exists
if [ ! -f "$BARCODE_FILE" ]; then
    echo "Error: Barcode file not found at: $BARCODE_FILE"
    exit 1
fi

# Create the destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

# Create a log file with a timestamp
LOG_DIR="${DEST_DIR}/logs"
mkdir -p "$LOG_DIR"
LOG_FILE="${LOG_DIR}/prep_fastq_$(date +%Y%m%d_%H%M%S).log"

# Redirect stdout and stderr to a log file and the console
exec &> >(tee -a "$LOG_FILE")

echo "--- Script Execution Started: $(date) ---"
echo "Source Directory: $SOURCE_DIR"
echo "Destination Directory: $DEST_DIR"
echo "Barcode Mapping File: $BARCODE_FILE"
echo "-------------------------------------------------"

# Read the mapping file line by line
while IFS=':' read -r filename barcode || [[ -n "$filename" ]]; do
    # Skip empty lines or lines where one part is missing
    if [ -z "$filename" ] || [ -z "$barcode" ]; then
        continue
    fi

    echo "Processing barcode: $barcode -> $filename.fastq.gz"
    
    BARCODE_SUBDIR="${SOURCE_DIR}/${barcode}"
    OUTPUT_FQ="${DEST_DIR}/${filename}.fastq.gz"

    # Check if the source barcode directory exists
    if [ ! -d "$BARCODE_SUBDIR" ]; then
        echo "WARNING: Directory does not exist, skipping: $BARCODE_SUBDIR"
        continue
    fi
    
    # Check if there are any gzipped fastq files to process
    if ! ls "${BARCODE_SUBDIR}"/*.fastq.gz 1> /dev/null 2>&1; then
        echo "WARNING: No *.fastq.gz files found in $BARCODE_SUBDIR, skipping."
        continue
    fi

    # Concatenate directly to the destination file
    cat "${BARCODE_SUBDIR}"/*.fastq.gz > "$OUTPUT_FQ"

    echo "Successfully created: $OUTPUT_FQ"
    echo "Finished processing barcode: $barcode"
    echo ""

done < "$BARCODE_FILE"

echo "-------------------------------------------------"
echo "--- Script Execution Finished: $(date) ---"
echo "Log file saved to: $LOG_FILE"