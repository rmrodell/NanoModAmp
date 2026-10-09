#!/usr/bin/bash
#SBATCH --job-name=20250503_endoBID_%j
#SBATCH --output=20250503_endoBID_%j.out
#SBATCH --error=20250503_endoBID_%j.err
#SBATCH --time=3:00:00
#SBATCH -p normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#SBATCH --mail-user=rodell@stanford.edu
#SBATCH --mail-type=BEGIN,END,FAIL

ml R


# fasta prep
bash /home/users/rodell/BIDamplicon/pipeline/fastaprep.sh /oak/stanford/groups/nicolemm/rodell/SHAPE-BID/20250418/20250418_1351_P2S-01916-A_PBC50989_3061744a/full_basecalling/pass /scratch/users/rodell/20250418_endoBID/20250503_70/fastq /scratch/users/rodell/20250418_endoBID/fastq/barcodes.txt
wait


# map to 4 psi, no deduplication
bash /home/users/rodell/BIDamplicon/pipeline/mapping_4psi.sh /scratch/users/rodell/20250418_endoBID/20250503_70/fastq /home/users/rodell/fastas/4psi.fasta /scratch/users/rodell/20250418_endoBID/20250503_70/4psi

# count 4psi
bash /home/users/rodell/BIDamplicon/pipeline/DelDetect.sh /scratch/users/rodell/20250418_endoBID/20250503_70/4psi/minimap2/sorted /scratch/users/rodell/20250418_endoBID/20250503_70/4psi/counts /home/users/rodell/fastas/4psi.fasta /home/users/rodell/4psi.bed

# deduplication and mapping
INPUT_DIR=/scratch/users/rodell/20250418_endoBID/20250503_70/fastq
REFERENCE_FASTA=/home/users/rodell/fastas/endoBIDamp_20250404.fasta
BASE_DIR=/scratch/users/rodell/20250418_endoBID/20250503_70/deduplication
THREADS=15

bash /home/users/rodell/BIDamplicon/pipeline/dedup/dedup_map_folder.sh "$INPUT_DIR" "$REFERENCE_FASTA" "$BASE_DIR" "$THREADS"

# count the deduplication
bash /home/users/rodell/BIDamplicon/pipeline/DelDetect.sh /scratch/users/rodell/20250418_endoBID/20250503_70/deduplication/minimap2/sorted /scratch/users/rodell/20250418_endoBID/20250503_70/deduplication/counts /home/users/rodell/fastas/endoBIDamp_20250404.fasta /scratch/users/rodell/20250404_endoBID/counts/20250404_endoBIDamplicon.bed
