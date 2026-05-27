#!/bin/bash
#SBATCH --job-name=fetch_hmm
#SBATCH --partition=express
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=fetch_hmm.out
#SBATCH --error=fetch_hmm.err


module load conda
conda activate hmmer_3.4

python hmm_fetch_1.py
#hmmfetch /home/lperez/TFG/4_alignment/4_1_run_hmmpress/data/Pfam-A.hmm 7tm_1 > /home/lperez/TFG/4_alignment/4_2_hmmfetch/results/PF00001.hmm

