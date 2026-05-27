# Postprocess Stockholm alignment files by reconstructing fragmented sequences and reformatting identifiers using UniProt information

import os
import csv

##### INPUT AND OUTPUT DIRECTORIES #####
# Input directory containing Stockholm alignment files from hmmalign
input_dir = "../4_1_Alignment_server/C_Hmmalign/results_sto"

# UniProt reference file for ID mapping (Entry -> Entry Name)
uniprot_tsv = "../../1_Table_uniprot_data/uniprotkb_organism_id_9606_AND_reviewed_2026_02_04.tsv"

# Output directory for reformatted alignment files
output_dir = "./results"
os.makedirs(output_dir, exist_ok=True)

ID_WIDTH = 30  # Fixed width so sequences start in the same column

##### LOAD UNIPROT MAPPING #####
# Create dictionary mapping UniProt Entry (accession) to Entry Name
uniprot_mapping = {}

with open(uniprot_tsv, "r") as f:
    reader = csv.DictReader(f, delimiter="\t")  # Read TSV using column names
    for row in reader:
        uniprot_mapping[row["Entry"]] = row["Entry Name"]  # Save mapping: Entry -> Entry Name

##### PROCESS ALL STOCKHOLM FILES #####
# Iterate over all .sto files in the input directory
for filename in os.listdir(input_dir):

    input_path = os.path.join(input_dir, filename)
    base_name = os.path.splitext(filename)[0]
    output_path = os.path.join(output_dir, base_name + ".txt")

    sequences = {}  # Dictionary to accumulate full aligned sequences per ID

    # Read Stockholm file and reconstruct fragmented sequences
    with open(input_path, "r") as f:
        for line in f:

            # Skip comment lines, empty lines, and end markers
            if line.lstrip().startswith("#") or line.strip() == "" or line.strip() == "//":
                continue

            # Split line into ID and aligned fragment
            parts = line.rstrip("\n").split(None, 1)
            if len(parts) < 2:
                continue

            raw_id, aligned_fragment = parts
            id_fields = raw_id.split("|")

            # Extract UniProt Entry and region coordinates
            if len(id_fields) >= 3:
                entry = id_fields[0]
                region = id_fields[-1]

                # Map Entry to Entry Name using UniProt reference
                entry_name = uniprot_mapping.get(entry, entry)
                new_id = f"{entry_name}/{region}"

                # Concatenate fragments if ID already exists (multi-line sequences)
                if new_id in sequences:
                    sequences[new_id] += aligned_fragment
                else:
                    sequences[new_id] = aligned_fragment

            else:
                continue

    # Write reformatted alignment file with one line per sequence
    with open(output_path, "w") as f:
        for new_id, full_sequence in sequences.items():
            f.write(f"{new_id.ljust(ID_WIDTH)} {full_sequence}\n")

print("All alignment files postprocessed successfully.")