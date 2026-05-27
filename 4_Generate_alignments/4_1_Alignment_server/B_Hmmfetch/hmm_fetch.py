# Extract individual HMM profiles from Pfam database using hmmfetch (HMMER suite)

import os
from subprocess import Popen


##### INPUT FILES #####
data_path = "/home/lperez/TFG/4_alignment/4_2_hmmfetch/data/pfam_code_short_name.tsv"
pfam_db = "/home/lperez/TFG/4_alignment/4_1_run_hmmpress/data/Pfam-A.hmm"
results_dir = "/home/lperez/TFG/4_alignment/4_2_hmmfetch/results"

os.makedirs(results_dir, exist_ok=True)


##### EXTRACT HMM PROFILES #####
# Open input file and iterate over Pfam codes
with open(data_path) as file:
    next(file)  # Skip header line

    for line in file:
        # Split line into Pfam code and short name
        pfam_code, short_name = line.rstrip("\n").split("\t")

        # Run hmmfetch to extract HMM profile for this domain
        with Popen(
            f"hmmfetch {pfam_db} {short_name} > {results_dir}/{pfam_code}.hmm",
            shell=True,
        ) as process:
            process.communicate()

print("HMM extraction completed successfully.")