#!/bin/bash
#SBATCH --job-name=align_hmm
#SBATCH --partition=express
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=align_hmm.out
#SBATCH --error=align_hmm.err


module load conda
conda activate hmmer_3.4

python 3_hmmalign_1.py
#hmmalign /home/lperez/TFG/4_alignment/4_2_hmmfetch/results/PF00001.hmm /home/lperez/TFG/4_alignment/4_3_hmmalign/domains_fasta/PF00001.fasta > /home/lperez/TFG/4_alignment/4_3_hmmalign/results_sto/PF00001_human_hmmalign.sto