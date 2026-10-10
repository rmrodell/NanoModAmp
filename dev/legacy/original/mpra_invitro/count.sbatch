#!/usr/bin/bash
#SBATCH --job-name=count_ivBID_%j
#SBATCH --output=count_ivBID_%j.out
#SBATCH --error=count_ivBID_%j.err
#SBATCH --time=1:00:00
#SBATCH -p normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --mail-user=rodell@stanford.edu
#SBATCH --mail-type=BEGIN,END,FAIL


bash $HOME/BIDamplicon/pipeline/BIDdetect.sh \
  --bam_dir /scratch/users/rodell/20251106_pool1/invitro/bam_links \
  --output_dir /scratch/users/rodell/20251106_pool1/invitro/counts_full \
  --ref_fasta /home/groups/nicolemm/rodell/pool1/pool1_cleaned_noadapters.fasta \
  --bed_file /home/groups/nicolemm/rodell/pool1/pool1_cleaned_full_noadapters.bed \
  --col_names 'celltype_vector_treat_rep'