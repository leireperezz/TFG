# Filter REC and SEED alignments to keep only shared protein sequences for fair comparison

import os

##### PATHS #####

# Directory containing REC alignments (postprocessed output)
rec_dir = "../4_Generate_alignments/4_3_Postprocess_alignments/results"

# Directory containing SEED alignments
seed_dir = "./data/Pfam_SEED_HUMAN_2024_03_20"

# Output directory for filtered REC alignments
output_rec_dir = "./results/rec_and_seed/REC_filtered"

# Output directory for filtered SEED alignments
output_seed_dir = "./results/rec_and_seed/SEED_filtered"

# Create output directories if they do not exist
os.makedirs(output_rec_dir, exist_ok=True)
os.makedirs(output_seed_dir, exist_ok=True)

##### CREATE PFAM TO FILE MAPPING #####

# Dictionary: PF00001 -> corresponding REC file name
rec_dict = {}
for f in os.listdir(rec_dir):
    if f.endswith(".txt"):
        pfam = f.split("_")[0]  # Extract PFAM code (e.g., PF00001)
        rec_dict[pfam] = f

# Dictionary: PF00001 -> corresponding SEED file name
seed_dict = {}
for f in os.listdir(seed_dir):
    if f.endswith(".txt"):
        pfam = f.split("_")[0]  # Extract PFAM code
        seed_dict[pfam] = f

# Identify PFAM families present in both datasets
common_pfams = set(rec_dict.keys()) & set(seed_dict.keys())

##### INITIALIZE COUNTER FOR EMPTY PFAMS #####

empty_count = 0  # Counts PFAMs with no shared sequences

##### PROCESS EACH PFAM FAMILY #####

for pfam in common_pfams:

    # Get corresponding file names
    rec_file = rec_dict[pfam]
    seed_file = seed_dict[pfam]

    # Build full file paths
    rec_path = os.path.join(rec_dir, rec_file)
    seed_path = os.path.join(seed_dir, seed_file)

    # Read REC alignment (ignore empty lines)
    with open(rec_path, "r", encoding="utf-8") as f:
        rec_lines = [line.rstrip("\n") for line in f if line.strip()]

    # Read SEED alignment
    with open(seed_path, "r", encoding="utf-8") as f:
        seed_lines = [line.rstrip("\n") for line in f if line.strip()]

    ##### EXTRACT PROTEIN IDENTIFIERS #####
    # Only the part before "/" is used (e.g., TRFR_HUMAN from TRFR_HUMAN/42-320)

    rec_codes = {line.split()[0].split("/")[0] for line in rec_lines}
    seed_codes = {line.split()[0].split("/")[0] for line in seed_lines}

    # Identify shared protein identifiers
    common_codes = rec_codes & seed_codes

    ##### FILTER ALIGNMENTS #####
    # Keep all lines whose protein code exists in both alignments

    rec_filtered = [
        line for line in rec_lines
        if line.split()[0].split("/")[0] in common_codes
    ]

    seed_filtered = [
        line for line in seed_lines
        if line.split()[0].split("/")[0] in common_codes
    ]

    ##### HANDLE EMPTY RESULTS #####
    # If no shared sequences remain, skip this PFAM

    if len(rec_filtered) == 0 or len(seed_filtered) == 0:
        empty_count += 1
        continue

    ##### SAVE FILTERED ALIGNMENTS #####

    # Save filtered REC file
    with open(os.path.join(output_rec_dir, rec_file), "w", encoding="utf-8") as f:
        f.write("\n".join(rec_filtered) + "\n")

    # Save filtered SEED file
    with open(os.path.join(output_seed_dir, seed_file), "w", encoding="utf-8") as f:
        f.write("\n".join(seed_filtered) + "\n")

##### FINAL OUTPUT #####

print(f"Number of empty PFAM files: {empty_count}")
print("Filtering completed successfully.")