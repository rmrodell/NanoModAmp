#!/bin/bash
set -e
set -u

ml biology
ml py-cutadapt/1.18_py36
ml samtools


if [ $# -lt 3 ]; then
    echo "Usage: $0 <input_fastq> <reference_fasta> <base_directory> [threads]"
    echo "Example: $0 /path/to/input.fq /path/to/reference.fasta /scratch/users/username/output_dir 7"
    echo "Note: If threads not specified, default value of 1 will be used"
    exit 1
fi

echo "--------------------------------------"

# Assign command line arguments to variables
INPUT_FILE=$1
REFERENCE_FASTA=$2
BASE_DIR=$3
THREADS=${4:-1}  # Use 1 threads if not specified


# Verify input files exist
if [ ! -f "${INPUT_FILE}" ]; then
    echo "Error: Input file ${INPUT_FILE} not found"
    exit 1
fi

if [ ! -f "${REFERENCE_FASTA}" ]; then
    echo "Error: Reference fasta ${REFERENCE_FASTA} not found"
    exit 1
fi

# Extract sample name from input file
SAMPLE_NAME=$(basename ${INPUT_FILE} .fq)
echo "Extracted sample name: ${SAMPLE_NAME}"

# Create directory structure
mkdir -p ${BASE_DIR}/{logs,trimmed1,UMI_extracted,minimap2,minimap2/sorted,minimap2/dedup}

echo "Processing ${SAMPLE_NAME}"
echo "Input file: ${INPUT_FILE}"
echo "Reference fasta: ${REFERENCE_FASTA}"
echo "Base directory: ${BASE_DIR}"




# Trim Adapters
trim_log=${BASE_DIR}/logs/trimmed1.log

echo "trimming Nanopore adapters"
cutadapt -m 125 -O 15 --cores ${THREADS} -a 'GAAGATAGAGCGACAGGCAAGT' -o ${BASE_DIR}/trimmed1/${SAMPLE_NAME}.fq ${INPUT_FILE} > $trim_log


#Extract UMI

r1=${BASE_DIR}/trimmed1/${SAMPLE_NAME}.fq
r1_umi=${BASE_DIR}/UMI_extracted/${SAMPLE_NAME}_UMI.fq
extract_log=${BASE_DIR}/logs/UMI_extract.log

echo "extracting UMIs"
umi_tools extract -I $r1 --extract-method=string --bc-pattern=NNNNNNNNNN --3prime -L $extract_log > $r1_umi


# Mapping
minimap2_log=${BASE_DIR}/logs/UMI_minimap2.log

echo "mapping with minimap2"
# Map using minimap2
/oak/stanford/groups/nicolemm/rodell/minimap2/minimap2 -a ${REFERENCE_FASTA} ${BASE_DIR}/UMI_extracted/${SAMPLE_NAME}_UMI.fq -k5 -t ${THREADS} > ${BASE_DIR}/minimap2/${SAMPLE_NAME}_UMI.sam 2>> ${minimap2_log}


# Convert to BAM file, sort the file, and index it
echo "Converting and sorting"
samtools view -b -o ${BASE_DIR}/minimap2/${SAMPLE_NAME}_UMI.bam ${BASE_DIR}/minimap2/${SAMPLE_NAME}_UMI.sam
samtools sort -@ ${THREADS} -O bam -o ${BASE_DIR}/minimap2/sorted/${SAMPLE_NAME}_UMI_sort.bam ${BASE_DIR}/minimap2/${SAMPLE_NAME}_UMI.bam
(cd ${BASE_DIR}/minimap2/sorted/ && samtools index ${SAMPLE_NAME}_UMI_sort.bam)


# deduplicate
echo "deduplicating UMIs"
umi_tools dedup --method directional -I ${BASE_DIR}/minimap2/sorted/${SAMPLE_NAME}_UMI_sort.bam -S ${BASE_DIR}/minimap2/dedup/${SAMPLE_NAME}_UMI_dedup.bam -L ${BASE_DIR}/logs/UMI_dedup.log

samtools index ${BASE_DIR}/minimap2/dedup/${SAMPLE_NAME}_UMI_dedup.bam

echo "summarizing read counts"
# Create summary file
SUMMARY_FILE=${BASE_DIR}/logs/${SAMPLE_NAME}_summary.txt
echo "----------------------------------------" >> ${SUMMARY_FILE}
echo "Processing Summary for ${SAMPLE_NAME}" > ${SUMMARY_FILE}
echo "----------------------------------------" >> ${SUMMARY_FILE}

# Count initial fastq reads
INITIAL_READS=$(echo $(cat ${INPUT_FILE} | wc -l)/4 | bc)
echo "Initial reads: ${INITIAL_READS}" >> ${SUMMARY_FILE}

# Count trimmed1 reads
TRIMMED1_READS=$(echo $(cat ${BASE_DIR}/trimmed1/${SAMPLE_NAME}.fq | wc -l)/4 | bc)
echo "After first trimming: ${TRIMMED1_READS}" >> ${SUMMARY_FILE}

# Count mapped reads (from BAM)
MAPPED_READS=$(samtools view -c ${BASE_DIR}/minimap2/sorted/${SAMPLE_NAME}_UMI_sort.bam)
echo "Mapped reads: ${MAPPED_READS}" >> ${SUMMARY_FILE}

# Count deduplicated reads
DEDUP_READS=$(samtools view -c ${BASE_DIR}/minimap2/dedup/${SAMPLE_NAME}_UMI_dedup.bam)
echo "After deduplication: ${DEDUP_READS}" >> ${SUMMARY_FILE}

# Calculate percentages

echo "Percentages relative to initial reads:" >> ${SUMMARY_FILE}
echo "First trimming: $(echo "scale=4; ${TRIMMED1_READS}/${INITIAL_READS}*100" | bc)%" >> ${SUMMARY_FILE}
echo "Mapped: $(echo "scale=4; ${MAPPED_READS}/${INITIAL_READS}*100" | bc)%" >> ${SUMMARY_FILE}
echo "After deduplication: $(echo "scale=4; ${DEDUP_READS}/${INITIAL_READS}*100" | bc)%" >> ${SUMMARY_FILE}

echo -e "\nPercentages relative to previous step:" >> ${SUMMARY_FILE}
echo "First trimming: $(echo "scale=2; ${TRIMMED1_READS}/${INITIAL_READS}*100" | bc)%" >> ${SUMMARY_FILE}
echo "Mapped: $(echo "scale=2; ${MAPPED_READS}/${TRIMMED1_READS}*100" | bc)%" >> ${SUMMARY_FILE}
echo "After deduplication: $(echo "scale=2; ${DEDUP_READS}/${MAPPED_READS}*100" | bc)%" >> ${SUMMARY_FILE}
