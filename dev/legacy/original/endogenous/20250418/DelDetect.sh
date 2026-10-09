#!/bin/bash
  
ml R

# Get the directory where DelDetect.sh is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# first argument is directory with find sorted bam files NO FINAL SLASH
# second argument is where outputs should go
# third argument is path to reference fasta file
# fourth argument is path to bed file with regions to examine


# Base path to the sorted & indexed bam files
BAM_PATH="$1"
DEST_PATH="$2"

# Name other input files
ref_fa="$3"
bed="$4"

# name where final counts will go
count_file="${DEST_PATH}/deldetect_counts.txt"

# Initialize the final output file with headers
echo -e "sample\tchr\tpos\tgene\ttotalReads\tA.count\tC.count\tG.count\tT.count\tDeletion.count\tInsertion.count\tref\tkmer\tstrand" > "$count_file"

(cd $BAM_PATH
# Loop through each .bam file in the input directory
for bam_file in "$BAM_PATH"/*.bam; do

    sample=$(basename "${bam_file}" _sort.bam)  # Get the base name without the extension

    echo "Processing sample: ${sample}"

    # extract counts using bam_counts.R
    Rscript "${SCRIPT_DIR}/bam_counts.R" --bedFile $bed --bamFile $bam_file --referenceFasta $ref_fa --outputFile ${sample}_counts.txt

    # Process the resultant file and append it to the final output
    result_file="${sample}_counts.txt"

    # Append the rest of the lines (skip header)
    awk -v sample="$sample" 'NR > 1 {print sample "\t" $0}' "$result_file" >> "$count_file"
done
)

echo "All files counted. Final count file located at $count_file" 

anova_input="${DEST_PATH}/deldetect_factors.txt"

# split sample name into different columns
Rscript "${SCRIPT_DIR}/sample_name.R" $count_file $anova_input celltype_vector_rep_treat
