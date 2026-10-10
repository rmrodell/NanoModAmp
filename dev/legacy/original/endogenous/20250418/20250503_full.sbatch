#!/usr/bin/bash
#SBATCH --job-name=20250503_endoBID_%j
#SBATCH --output=20250503_endoBID_%j.out
#SBATCH --error=20250503_endoBID_%j.err
#SBATCH --time=6:00:00
#SBATCH -p normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=64GB
#SBATCH --mail-user=rodell@stanford.edu
#SBATCH --mail-type=BEGIN,END,FAIL

ml R

bash /home/users/rodell/BIDamplicon/pipeline/DelDetect.sh /scratch/users/rodell/20250418_endoBID/20250503_70/deduplication/minimap2/sorted /scratch/users/rodell/20250418_endoBID/20250503_70/deduplication/counts /home/users/rodell/fastas/endoBIDamp_20250404.fasta /scratch/users/rodell/20250418_endoBID/20250503_70/20250404_endoBIDamplicon_full.bed

