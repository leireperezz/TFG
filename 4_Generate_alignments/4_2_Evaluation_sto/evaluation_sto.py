# Evaluate alignment generation by comparing FASTA input sequences with Stockholm output files

import pandas as pd
import os

##### INPUT PATHS #####
# FASTA files from domain extraction step
fasta_path = "../../3_Extract_sequences/results/domains_fasta"

# Stockholm alignment files from hmmalign step
sto_path = "../4_1_Alignment_server/C_Hmmalign/results_sto"

##### COUNT SEQUENCES IN FASTA FILES #####
fasta_counts = {}

# Iterate over FASTA files and count sequences (headers starting with >)
for filename in os.listdir(fasta_path):
    if not filename.lower().endswith((".fasta", ".fa", ".faa")):
        continue
    
    pfam_id = os.path.splitext(filename)[0]
    file_path = os.path.join(fasta_path, filename)

    num_seqs = 0

    with open(file_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                num_seqs = num_seqs + 1

    fasta_counts[pfam_id] = num_seqs

print(f"FASTAs processed: {len(fasta_counts)}")

# Count Pfam domains with only 1 sequence
pfams_1_seq = [pfam for pfam, n in fasta_counts.items() if n == 1]
print(f"Pfams with 1 sequence: {len(pfams_1_seq)}")

##### COUNT LINES IN STOCKHOLM FILES #####
sto_line_counts = {}

# Iterate over Stockholm files and count total lines
for filename in os.listdir(sto_path):
    if filename.endswith(".sto"):
        
        pfam_id = filename.split("_")[0]
        file_path = os.path.join(sto_path, filename)

        with open(file_path, "r") as f:
            num_lines = sum(1 for _ in f)

        sto_line_counts[pfam_id] = num_lines

print(f"STOs processed: {len(sto_line_counts)}")

# Count empty Stockholm files (0 lines)
num_empty = sum(1 for v in sto_line_counts.values() if v == 0)
print(f"STOs with 0 lines: {num_empty}")

##### MERGE DATA #####
# Combine all Pfam IDs from both FASTA and STO files
all_pfams = set(fasta_counts.keys()) | set(sto_line_counts.keys())

data = []

# Create summary table with FASTA sequence counts and STO line counts
for pfam in sorted(all_pfams):
    data.append({
        "pfam_id": pfam,
        "fasta_seqs": fasta_counts.get(pfam, 0),
        "sto_lines": sto_line_counts.get(pfam, 0)
    })

df_summary = pd.DataFrame(data)

##### SAVE COMPLETE SUMMARY #####
output_dir = "./results"
os.makedirs(output_dir, exist_ok=True)

output_file = os.path.join(output_dir, "sto_fasta_summary.csv")

df_summary.to_csv(output_file, index=False)
print(f"Summary file saved at: {output_file}")

##### SAVE FILTERED SUMMARY (FASTA = 1 SEQUENCE) #####
df_only_1 = df_summary[df_summary["fasta_seqs"] == 1]
output_file_1 = os.path.join(output_dir, "sto_fasta_summary_only_1_seq.csv")
df_only_1.to_csv(output_file_1, index=False)
print(f"Filtered file (1 seq) saved at: {output_file_1}")

##### SAVE FILTERED SUMMARY (STO = 0 LINES) #####
df_sto_0 = df_summary[df_summary["sto_lines"] == 0]
output_file_0 = os.path.join(output_dir, "sto_fasta_summary_only_sto_0.csv")
df_sto_0.to_csv(output_file_0, index=False)
print(f"Filtered file (0 lines) saved at: {output_file_0}")