#!/usr/bin/bash
#SBATCH --job-name=count_icBID_%j
#SBATCH --output=count_icBID_%j.out
#SBATCH --error=count_icBID_%j.err
#SBATCH --time=4:00:00
#SBATCH -p normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --mail-user=rodell@stanford.edu
#SBATCH --mail-type=BEGIN,END,FAIL

bash $HOME/BIDamplicon/pipeline/BIDdetect.sh \
  --bam_dir /scratch/users/rodell/20251106_pool1/incell/prep/final_bam \
  --output_dir /scratch/users/rodell/20251106_pool1/incell/counts_full \
  --ref_fasta /home/groups/nicolemm/rodell/pool1/pool1_cleaned_noadapters.fasta \
  --bed_file /home/groups/nicolemm/rodell/pool1/pool1_cleaned_full_noadapters.bed
