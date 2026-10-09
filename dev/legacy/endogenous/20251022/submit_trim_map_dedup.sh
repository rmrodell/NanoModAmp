#!/bin/bash
#
# SLURM job array submitter for a general bioinformatics pipeline.
# This script takes all configuration from the command line, making it highly reusable.
#

# Stop the script if any command fails or if an unset variable is used
set -eu -o pipefail

################################################################################
### SCRIPT USAGE FUNCTION                                                  ###
################################################################################

# Function to display help message
usage() {
    echo "Usage: $0 -u <email> -s <script> -a <map> -i <input_dir> -o <out_dir> -r <ref_fa> [OPTIONS]..."
    echo ""
    echo "Submits a SLURM job array where each task processes one sample from the SAMPLE_MAP."
    echo ""
    echo "Required Arguments:"
    echo "  -u, --mail-user <email>      Email for SLURM notifications."
    echo "  -s, --script-path <path>     Path to the executable pipeline script (e.g., process_single_file.sh)."
    echo "  -a, --sample-map <path>      Path to the text file listing samples, one per line."
    echo "  -i, --input-dir <path>       Path to the directory containing the source FASTQ files."
    echo "  -o, --output-dir <path>      Top-level directory where all outputs will be created."
    echo "  -r, --ref_fa <path>          Reference fasta for mapping files."
    echo ""
    echo "Optional Job Resource Arguments (with defaults):"
    echo "  -c, --cpu <int>              CPUs per job task. (Default: 16)"
    echo "  -m, --mem <string>           Memory per job task (e.g., 32G). (Default: 32G)"
    echo "  -t, --time <string>          Time limit for each job task (e.g., 24:00:00). (Default: 24:00:00)"
    echo "  -p, --partition <string>     SLURM partition to use. (Default: normal)"
    echo "  -h, --help                   Display this help message and exit."
    echo ""
    echo "Example:"
    echo "  $0 \\"
    echo "      --mail-user me@stanford.edu \\"
    echo "      --script-path /path/to/process_single_file.sh \\"
    echo "      --sample-map /path/to/my_samples.txt \\"
    echo "      --input-dir /path/to/fastq_files/ \\"
    echo "      --output-dir /path/to/my_output_project/ \\"
    echo "      --cpu 8 --mem 16G"
}

################################################################################
### DEFAULT SETTINGS & ARGUMENT PARSING                                    ###
################################################################################

# --- Set Default Optional SLURM Job Resource Configuration ---
CPU_PER_JOB=16
MEM_PER_JOB="32G"
TIME_LIMIT="12:00:00"
PARTITION="normal"

# --- Initialize Required Argument Variables ---
MAIL_USER=""
PIPELINE_SCRIPT_PATH=""
SAMPLE_MAP_FILE=""
INPUT_SOURCE_DIR=""
TOP_LEVEL_OUTPUT_DIR=""

# --- Parse All Arguments ---
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -u|--mail-user)         MAIL_USER="$2"; shift 2 ;;
        -s|--script-path)       PIPELINE_SCRIPT_PATH="$2"; shift 2 ;;
        -a|--sample-map)        SAMPLE_MAP_FILE="$2"; shift 2 ;;
        -i|--input-dir)         INPUT_SOURCE_DIR="$2"; shift 2 ;;
        -o|--output-dir)        TOP_LEVEL_OUTPUT_DIR="$2"; shift 2 ;;
        -r|--ref_fa)            REF_FA="$2"; shift 2 ;;
        -c|--cpu)               CPU_PER_JOB="$2"; shift 2 ;;
        -m|--mem)               MEM_PER_JOB="$2"; shift 2 ;;
        -t|--time)              TIME_LIMIT="$2"; shift 2 ;;
        -p|--partition)         PARTITION="$2"; shift 2 ;;
        -h|--help)              usage; exit 0 ;;
        *)                      echo "Error: Unknown option '$1'"; usage; exit 1 ;;
    esac
done

# --- Final Input Validation ---
# Check for missing required arguments
missing_args=()
if [ -z "$MAIL_USER" ]; then missing_args+=("--mail-user"); fi
if [ -z "$PIPELINE_SCRIPT_PATH" ]; then missing_args+=("--script-path"); fi
if [ -z "$SAMPLE_MAP_FILE" ]; then missing_args+=("--sample-map"); fi
if [ -z "$INPUT_SOURCE_DIR" ]; then missing_args+=("--input-dir"); fi
if [ -z "$TOP_LEVEL_OUTPUT_DIR" ]; then missing_args+=("--output-dir"); fi
if [ -z "$REF_FA" ]; then missing_args+=("--ref_fa"); fi

if [ ${#missing_args[@]} -gt 0 ]; then
    echo "Error: The following required arguments were not provided:"
    for arg in "${missing_args[@]}"; do
        echo "  $arg"
    done
    usage
    exit 1
fi

# Check that paths/files actually exist
if [ ! -x "$PIPELINE_SCRIPT_PATH" ]; then echo "Error: Pipeline script '$PIPELINE_SCRIPT_PATH' not found or is not executable."; exit 1; fi
if [ ! -f "$SAMPLE_MAP_FILE" ]; then echo "Error: Sample map file '$SAMPLE_MAP_FILE' not found."; exit 1; fi
if [ ! -d "$INPUT_SOURCE_DIR" ]; then echo "Error: FASTQ source directory '$INPUT_SOURCE_DIR' not found."; exit 1; fi
if [ ! -r "$REF_FA" ]; then echo "Error: Reference fasta '$REF_FA' not found."; exit 1; fi


################################################################################
### BATCH SUBMISSION LOGIC                                                 ###
################################################################################

echo "--- Preparing SLURM Job Array Submission ---"
echo "  CPUs per Task:    $CPU_PER_JOB"
echo "  Memory per Task:  $MEM_PER_JOB"
echo "  Time Limit:       $TIME_LIMIT"
echo "  Partition:        $PARTITION"
echo "  Notify Email:     $MAIL_USER"
echo "  Pipeline Script:  $PIPELINE_SCRIPT_PATH"
echo "  Sample Map:       $SAMPLE_MAP_FILE"
echo "  FASTQ Source Dir: $INPUT_SOURCE_DIR"
echo "  Output Dir:       $TOP_LEVEL_OUTPUT_DIR"
echo "  Reference Fasta:  $REF_FA"
echo "---------------------------------------------"

# Define and create output directories
SLURM_LOGS_DIR="${TOP_LEVEL_OUTPUT_DIR}/slurm_logs"
mkdir -p "$SLURM_LOGS_DIR"

echo "SLURM output/error logs will be saved in: $SLURM_LOGS_DIR"

# Determine the number of jobs for the array
N_JOBS=$(wc -l < "$SAMPLE_MAP_FILE")
if [ "$N_JOBS" -eq 0 ]; then
    echo "Error: Sample map file '$SAMPLE_MAP_FILE' is empty."
    exit 1
fi
echo "Found $N_JOBS samples to process."

# Submit the array job
echo "Submitting array job..."
sbatch \
    --job-name="SHAPE_prep" \
    --array=1-"$N_JOBS" \
    --cpus-per-task="$CPU_PER_JOB" \
    --mem="$MEM_PER_JOB" \
    --time="$TIME_LIMIT" \
    --partition="$PARTITION" \
    --mail-type=BEGIN,END,FAIL \
    --mail-user="$MAIL_USER" \
    --output="${SLURM_LOGS_DIR}/job_%A_task_%a.out" \
    --error="${SLURM_LOGS_DIR}/job_%A_task_%a.err" \
    "$PIPELINE_SCRIPT_PATH" --map_file "$SAMPLE_MAP_FILE" --input_dir "$INPUT_SOURCE_DIR" --output_dir "$TOP_LEVEL_OUTPUT_DIR" --ref_fasta "$REF_FA"

echo "--- Job array submitted successfully! ---"
