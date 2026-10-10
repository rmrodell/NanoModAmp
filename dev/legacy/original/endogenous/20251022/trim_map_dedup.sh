#!/bin/bash
#
# A generalized pipeline for processing a SINGLE SHAPE-MaP sample.
# This script is designed to be called as a SLURM array task.
#

# --- 0. Argument Parsing with Flags ---

# Function to display usage information
usage() {
    echo "Usage: $0 -m <sample_map_file> -i <input_dir> -o <top_level_output_dir> -r <ref_fasta>"
    echo ""
    echo "Options:"
    echo "  -m, --map_file <file>        Path to the file mapping sample IDs to barcodes."
    echo "  -i, --input_dir <dir>        Directory containing the input fastq.gz files."
    echo "  -o, --output_dir <dir>       The root directory where all output will be stored."
    echo "  -r, --ref_fasta <file>       Reference FASTA file used for mapping."
    echo "  -h, --help                   Display this help message."
    exit 1
}

# Define short and long options
SHORT_OPTS="m:i:o:r:h"
LONG_OPTS="map_file:,input_dir:,output_dir:,ref_fasta:,help"

# Parse the options using getopt
PARSED_OPTS=$(getopt -o "$SHORT_OPTS" --long "$LONG_OPTS" -n "$0" -- "$@")

# Check if getopt had a problem
if [ $? -ne 0 ]; then
    echo "Error parsing options." >&2
    usage
fi

# Replace the script's positional parameters with the parsed options
eval set -- "$PARSED_OPTS"

# Loop through the parsed options to assign variables
while true; do
    case "$1" in
        -m|--map_file)
            SAMPLE_MAP_FILE="$2"
            shift 2
            ;;
        -i|--input_dir)
            FQ_SOURCE_DIR="$2"
            shift 2
            ;;
        -o|--output_dir)
            TOP_LEVEL_OUTPUT_DIR="$2"
            shift 2
            ;;
        -r|--ref_fasta)
            REF_FA="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Internal error!"
            exit 1
            ;;
    esac
done

# Check if all required arguments were provided
if [ -z "$SAMPLE_MAP_FILE" ] || [ -z "$FQ_SOURCE_DIR" ] || [ -z "$TOP_LEVEL_OUTPUT_DIR" ] || [ -z "$REF_FA" ]; then
    echo "Error: Missing one or more required arguments."
    usage
fi

# Check if input files/directories exist
if [ ! -f "$SAMPLE_MAP_FILE" ]; then
    echo "Error: Sample map file not found at: $SAMPLE_MAP_FILE"
    exit 1
fi
if [ ! -f "$REF_FA" ]; then
    echo "Error: Reference FASTA file not found at: $REF_FA"
    exit 1
fi
if [ ! -d "$FQ_SOURCE_DIR" ]; then
    echo "Error: FASTQ source directory not found at: $FQ_SOURCE_DIR"
    exit 1
fi


# Echo the assigned arguments to the log for traceability
echo "--- Script Arguments Received ---"
echo "Sample Map File:          ${SAMPLE_MAP_FILE}"
echo "FASTQ Source Directory:   ${FQ_SOURCE_DIR}"
echo "Top-Level Output Dir:     ${TOP_LEVEL_OUTPUT_DIR}"
echo "Reference FASTA:          ${REF_FA}"
echo "---------------------------------"
echo "" # Add a blank line for readability

# --- 1. SLURM Array Task Setup ---
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "Error: This script must be run as a SLURM array job."
    exit 1
fi

# Read the line from the sample map corresponding to the task ID
TASK_INFO=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_MAP_FILE")
if [ -z "$TASK_INFO" ]; then
    echo "Error: No data found in '$SAMPLE_MAP_FILE' for task ID $SLURM_ARRAY_TASK_ID."
    exit 1
fi

# Parse the sample name and barcode from the line
IFS=':' read -r SAMPLE_ID BARCODE <<< "$TASK_INFO"
SAMPLE_ID=$(echo "$SAMPLE_ID" | xargs)
BARCODE=$(echo "$BARCODE" | xargs)

# Find the specific input FASTQ file
RAW_FASTQ=$(find "$FQ_SOURCE_DIR" -name "${SAMPLE_ID}.fastq.gz" | head -n 1)
if [ -z "$RAW_FASTQ" ]; then
    echo "Error: Could not find FASTQ file for sample '$SAMPLE_ID' in '$FQ_SOURCE_DIR'."
    exit 1
fi

# Define output directories based on parsed info
PROJECT_DIR="${TOP_LEVEL_OUTPUT_DIR}/${SAMPLE_ID}"
FINAL_COMMON_DIR="${TOP_LEVEL_OUTPUT_DIR}/deduplicated_bam"

# Load modules
ml biology py-cutadapt/1.18_py36 samtools
ml python/3.6.1

# Stop the script if any command fails
set -e
set -o pipefail

################################################################################
###                  USER DEFINED VARIABLES                                  ###
################################################################################
SENSE_ADAPTER_5PRIME="TTTCTGTTGGTGCTGATATTGCG"
SENSE_ADAPTER_3PRIME="GAAGATAGAGCGACAGGCAAGT"
ANTISENSE_ADAPTER_5PRIME="ACTTGCCTGTCGCTCTATCTTC"
ANTISENSE_ADAPTER_3PRIME="CGCAATATCAGCACCAACAGAAA"
UMI_PATTERN="NNNNNNNNNN"

# SLURM-aware thread count
if [ -n "$SLURM_CPUS_PER_TASK" ]; then THREADS="$SLURM_CPUS_PER_TASK"; else THREADS=1; fi

################################################################################
###             HELPER FUNCTIONS AND DIRECTORY SETUP                         ###
################################################################################
log_message() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [Job:${SLURM_ARRAY_JOB_ID} Task:${SLURM_ARRAY_TASK_ID}] [$SAMPLE_ID] $1"; }
format_duration() { local s=$1; echo "$((s/60))m $((s%60))s"; }

track_metrics() {
    local step_name=$1; local file_path=$2
    if [ ! -s "$file_path" ]; then
        echo -e "${step_name}\t${SAMPLE_ID}\t$(basename "$file_path")\t0\t0B" >> "$METRICS_FILE"; return
    fi
    local size=$(du -h "$file_path" | cut -f1); local count
    if [[ "$file_path" == *.fastq || "$file_path" == *.fq ]]; then
        count=$(echo "$(wc -l < "$file_path") / 4" | bc)
    elif [[ "$file_path" == *.bam || "$file_path" == *.sam ]]; then
        count=$(samtools view -@ "$THREADS" -c "$file_path")
    else
        count="N/A"
    fi
    echo -e "${step_name}\t${SAMPLE_ID}\t$(basename "$file_path")\t${count}\t${size}" >> "$METRICS_FILE"
}

# --- CORRECTED & SIMPLIFIED DIRECTORY DEFINITIONS ---
TMP_DIR="${PROJECT_DIR}/tmp"
FINAL_DIR="${PROJECT_DIR}/final"
LOGS_DIR="${PROJECT_DIR}/logs"
REPORTS_DIR="${PROJECT_DIR}/reports"

mkdir -p "$TMP_DIR" "$FINAL_DIR" "$LOGS_DIR" "$REPORTS_DIR" "$FINAL_COMMON_DIR"

METRICS_FILE="${REPORTS_DIR}/${SAMPLE_ID}.run_metrics.tsv"; PIPELINE_LOG="${LOGS_DIR}/${SAMPLE_ID}.pipeline_run.log"
exec > >(tee -a "${PIPELINE_LOG}") 2>&1
echo -e "Step\tSample\tFile\tRecord_Count\tFile_Size" > "$METRICS_FILE"

################################################################################
### PIPELINE STEPS                                                           ###
################################################################################

log_message "--- Starting pipeline ---"
log_message "--- Using ${THREADS} threads for parallel tasks ---"
log_message "--- Input FASTQ: ${RAW_FASTQ} ---"
log_message "--- Project Dir: ${PROJECT_DIR} ---"

# --- Define intermediate file paths ---
SENSE_TRIM1_RAW="${TMP_DIR}/${SAMPLE_ID}_sense_trim1_raw.fastq"
ANTISENSE_TRIM1_RAW="${TMP_DIR}/${SAMPLE_ID}_antisense_trim1_raw.fastq"
ANTISENSE_TRIM1_RC="${TMP_DIR}/${SAMPLE_ID}_antisense_trim1_rc.fastq"
COMBINED_TRIM1="${TMP_DIR}/${SAMPLE_ID}_all_sense_trim1.fastq"
COMBINED_UMI="${TMP_DIR}/${SAMPLE_ID}_all_sense_umi.fastq"
MAPPED_SAM="${TMP_DIR}/${SAMPLE_ID}_mapped.sam"
MAPPED_BAM="${TMP_DIR}/${SAMPLE_ID}_mapped_sorted.bam"
DEDUP_BAM="${FINAL_DIR}/${SAMPLE_ID}_deduplicated.bam"
MINIMAP2_LOG="${LOGS_DIR}/${SAMPLE_ID}_minimap2_mapping.log"
DEDUP_LOG="${LOGS_DIR}/${SAMPLE_ID}_umi_tools_dedup.log"

# --- 1. Setup Read Tracking ---
log_message "Step 1: Start Read Tracking..."
start_time=$(date +%s)
track_metrics "1_Initial_FASTQ" "$RAW_FASTQ"
end_time=$(date +%s); log_message "Step 1 finished. Duration: $(format_duration $((end_time - start_time)))"

# --- 2. Filter and Trim FASTQ streams ---
log_message "Step 2: Filtering, trimming, and sanitizing streams..."
start_time=$(date +%s)

# Pass 2a: Filter Sense reads
log_message "  - Pass 2a: Filtering for Sense reads..."
cutadapt \
    -m 125 \
    -O 15 \
    --cores="$THREADS" \
    --discard-untrimmed \
    -g "$SENSE_ADAPTER_5PRIME...$SENSE_ADAPTER_3PRIME" \
    -o "$SENSE_TRIM1_RAW" \
    "$RAW_FASTQ" \
    > "${LOGS_DIR}/${SAMPLE_ID}_cutadapt_sense_filter.log"
track_metrics "2a_Sense_Trim1_Raw" "$SENSE_TRIM1_RAW" 

# Pass 2b: Filter Antisense reads
log_message "  - Pass 2b: Filtering for Antisense reads..."
cutadapt \
    -m 125 \
    -O 15 \
    --cores="$THREADS" \
    --discard-untrimmed \
    -g "$ANTISENSE_ADAPTER_5PRIME...$ANTISENSE_ADAPTER_3PRIME" \
    -o "$ANTISENSE_TRIM1_RAW" \
    "$RAW_FASTQ" \
    > "${LOGS_DIR}/${SAMPLE_ID}_cutadapt_antisense_filter.log"
track_metrics "2b_Antisense_Trim1_Raw" "$ANTISENSE_TRIM1_RAW" 

end_time=$(date +%s); log_message "Step 2 finished. Duration: $(format_duration $((end_time - start_time)))"


# --- 3. Unify Read Orientation ---
log_message "Step 3: Unifying read orientation..."
start_time=$(date +%s)
log_message "  - Reverse-complementing antisense reads..."

$HOME/seqtk/seqtk seq -r "$ANTISENSE_TRIM1_RAW" > "$ANTISENSE_TRIM1_RC"
track_metrics "3a_Antisense_RC" "$ANTISENSE_TRIM1_RC"
log_message "  - Combining all reads into a single sense-oriented stream..."

# Combine the sense file and the RC antisense file
cat "$SENSE_TRIM1_RAW" "$ANTISENSE_TRIM1_RC" > "$COMBINED_TRIM1"
track_metrics "3b_Combined_Trim1" "$COMBINED_TRIM1"
end_time=$(date +%s); log_message "Step 3 finished. Duration: $(format_duration $((end_time - start_time)))"

# --- 4. Extract UMIs from Combined File ---
log_message "Step 4: Extracting UMIs from 3' end of all reads..."
start_time=$(date +%s)
umi_tools extract -I "$COMBINED_TRIM1" --extract-method=string --bc-pattern="$UMI_PATTERN" --3prime -L "${LOGS_DIR}/${SAMPLE_ID}_umi_tools.log" > "$COMBINED_UMI"
track_metrics "4_Combined_UMI" "$COMBINED_UMI"
end_time=$(date +%s); log_message "Step 4 finished. Duration: $(format_duration $((end_time - start_time)))"

# --- 5. Map reads with Minimap2 ---
log_message "Step 5: Mapping reads with Minimap2..."
start_time=$(date +%s)

$OAK/rodell/minimap2/minimap2 \
    -a \
    "$REF_FA" \
    "$COMBINED_UMI" \
    -k5 -t "$THREADS" \
    > "$MAPPED_SAM"

track_metrics "5_Mapped_SAM" "$MAPPED_SAM"
end_time=$(date +%s); log_message "Step 5 finished. Duration: $(format_duration $((end_time - start_time)))"

# --- 6. Convert to Sorted BAM and Index ---
log_message "Step 6: Converting SAM to sorted BAM and indexing..."
start_time=$(date +%s)
samtools view -@ "$THREADS" -b "$MAPPED_SAM" | samtools sort -@ "$THREADS" -o "$MAPPED_BAM"
samtools index "$MAPPED_BAM"
track_metrics "6_Sorted_BAM" "$MAPPED_BAM"
end_time=$(date +%s); log_message "Step 6 finished. Duration: $(format_duration $((end_time - start_time)))"

# --- 7. Deduplicate Reads with UMI-tools ---
log_message "Step 7: Deduplicating reads with UMI-tools and indexing..."
start_time=$(date +%s)
umi_tools dedup \
    --method directional \
    -I "$MAPPED_BAM" \
    -S "$DEDUP_BAM" \
    -L "$DEDUP_LOG"
samtools index "$DEDUP_BAM"
track_metrics "7_Deduplicated_BAM" "$DEDUP_BAM"
end_time=$(date +%s); log_message "Step 7 finished. Duration: $(format_duration $((end_time - start_time)))"


# --- Step 8: Copy Final BAM to Common Output Directory ---
log_message "Step 8: Copying and indexing final BAM..."
start_time=$(date +%s)

DEST_BAM_PATH="${FINAL_COMMON_DIR}/${SAMPLE_ID}.bam"

# 1. Copy the file
cp "$DEDUP_BAM" "$DEST_BAM_PATH"
log_message "  - Copied to: ${DEST_BAM_PATH}"

# 2. Verify the copy with md5sum
log_message "  - Verifying file integrity with md5sum..."
source_md5=$(md5sum "$DEDUP_BAM" | awk '{print $1}')
dest_md5=$(md5sum "$DEST_BAM_PATH" | awk '{print $1}')

# 3. Compare checksums and act accordingly
if [ "$source_md5" == "$dest_md5" ]; then
    log_message "  - SUCCESS: Checksums match. File integrity confirmed."
else
    log_message "  - FATAL ERROR: Checksums DO NOT MATCH!"
    log_message "    - Source:      ${source_md5}"
    log_message "    - Destination: ${dest_md5}"
    log_message "  - Halting pipeline to prevent use of corrupted data."
    exit 1
fi

# 4. Index the bam file
samtools index "$DEST_BAM_PATH"

# --- Step 9: Clean up intermediate files ---
log_message "Step 9: Intermediate file cleanup is currently disabled."
# rm -rf "$TMP_DIR"

log_message "--- Pipeline for ${SAMPLE_ID} finished successfully! ---"

# --- Step 10: Display Final Metrics Report ---
log_message "Step 10: Displaying run metrics summary for ${SAMPLE_ID}:"
column -t -s $'\t' "$METRICS_FILE"