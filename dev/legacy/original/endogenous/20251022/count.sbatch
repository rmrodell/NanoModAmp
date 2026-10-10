#!/usr/bin/bash
#SBATCH --job-name=count_endoBID_%j
#SBATCH --output=count_endoBID_%j.out
#SBATCH --error=count_endoBID_%j.err
#SBATCH --time=6:00:00
#SBATCH -p normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --mail-user=rodell@stanford.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# bash $HOME/BIDamplicon/pipeline/BIDdetect.sh \
#  --bam_dir /scratch/users/rodell/20251022_endoBID/mapping/deduplicated_bam \
#  --output_dir /scratch/users/rodell/20251022_endoBID/counts_delpos \
#  --ref_fasta /scratch/users/rodell/20251022_endoBID/set3.fa \
#  --bed_file /scratch/users/rodell/20251022_endoBID/set3_delpos.bed

bash $HOME/BIDamplicon/pipeline/BIDdetect.sh \
  --bam_dir /scratch/users/rodell/20251022_endoBID/mapping/deduplicated_bam \
  --output_dir /scratch/users/rodell/20251022_endoBID/counts_full \
  --ref_fasta /scratch/users/rodell/20251022_endoBID/set3.fa \
  --bed_file /scratch/users/rodell/20251022_endoBID/set3_fulllength.bed