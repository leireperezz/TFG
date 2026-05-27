# Generate multiple sequence alignments using hmmalign (HMMER suite) for each Pfam domain

import os
from subprocess import Popen

##### INPUT FILES #####
fasta_dir = "/home/lperez/TFG/4_alignment/4_3_hmmalign/domains_fasta"
hmm_dir = "/home/lperez/TFG/4_alignment/4_2_hmmfetch/results"
output_dir = "/home/lperez/TFG/4_alignment/4_3_hmmalign/results_sto"

##### GENERATE ALIGNMENTS #####
# Iterate over all FASTA files in the input directory
for file in os.listdir(fasta_dir):
    if file.endswith(".fasta"):
        
        # Extract Pfam code from filename
        pfam_code = file.replace(".fasta", "")
        
        # Define file paths for HMM profile, FASTA sequences, and output alignment
        hmm_file = f"{hmm_dir}/{pfam_code}.hmm"
        fasta_file = f"{fasta_dir}/{file}"
        output_sto = f"{output_dir}/{pfam_code}_human_hmmalign.sto"

        # Run hmmalign to generate alignment in Stockholm format
        with Popen(
            f"hmmalign {hmm_file} {fasta_file} > {output_sto}",
            shell=True,
        ) as process:
            process.communicate()

print("Alignment generation completed successfully.")