#!/bin/bash
#SBATCH --job-name=pfam_press
#SBATCH --partition=express
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=pfam_press.out
#SBATCH --error=pfam_press.err

module load conda
conda activate hmmer_3.4

hmmpress ./data/Pfam-A.hmm  # run only once

# This command is needed to create indexed binary files (.h3f, .h3i, .h3m, .h3p):
# Quickly access individual profiles
# Use hmmfetch
# Run hmmsearch efficiently

