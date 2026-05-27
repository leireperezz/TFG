# Extract Pfam domain sequences from human proteome using coordinates obtained from InterPro API.
# Creates one FASTA file per Pfam domain and generates distribution statistics.

import pandas as pd  
import os            



##### INPUT FILES #####
# Read Pfam domain coordinates from previous step
pfam_file = "../2_Get_coordinates_pfam_domains/pfam_domains_all.tsv"
pfam_file = os.path.normpath(pfam_file)
df_pfam = pd.read_csv(pfam_file, sep="\t")

# Read UniProt protein sequences
uniprot_file = "../1_Table_uniprot_data/uniprotkb_organism_id_9606_AND_reviewed_2026_02_04.tsv"
uniprot_file = os.path.normpath(uniprot_file)
df_uniprot = pd.read_csv(uniprot_file, sep="\t")

# Rename column to use the same key in both tables
df_uniprot = df_uniprot.rename(columns={"Entry": "uniprot"})



##### RESULTS OUTPUT FOLDERS #####
# Get the directory where this script is located
script_dir = os.path.dirname(os.path.abspath(__file__))

# Create "results" folder in the script directory (if it does not exist)
results_dir = os.path.join(script_dir, "results")
os.makedirs(results_dir, exist_ok=True)

# Create "domains_fasta" subfolder inside results
domains_dir = os.path.join(results_dir, "domains_fasta")
os.makedirs(domains_dir, exist_ok=True)


##### DATA EXPLORATION #####
# Total number of unique UniProt proteins
total_proteins = df_uniprot["uniprot"].nunique()
print("Total unique UniProt proteins:", total_proteins)

# Number of proteins with at least one Pfam domain
proteins_with_pfam = df_pfam.loc[df_pfam["domain_id"].notna(), "uniprot"].nunique()
print("Proteins with at least one Pfam domain:", proteins_with_pfam)

# Number of rows without Pfam domain annotation
na_count = df_pfam["domain_id"].isna().sum()
print("Number of NA rows in domain_id:", na_count)

# Total number of distinct Pfam domain IDs
num_unique_domains = df_pfam["domain_id"].nunique()
print("Number of unique Pfam domain IDs:", num_unique_domains)



##### MERGE AND DOMAIN EXTRACTION #####
# Merge Pfam domain coordinates with UniProt protein sequences
df_merged = df_pfam.merge(
    df_uniprot[["uniprot", "Sequence"]],  # Keep only necessary columns
    on="uniprot",                         # Merge by UniProt ID
    how="left"                            # Keep all Pfam rows
)

# Convert start and end positions to nullable integers
df_merged["start"] = df_merged["start"].astype("Int64")
df_merged["end"] = df_merged["end"].astype("Int64")

# Extract domain subsequence using start/end coordinates (convert 1-based to 0-based indexing)
df_merged["domain_sequence"] = df_merged.apply(
    lambda row: (
        row["Sequence"][row["start"] - 1 : row["end"]]
        if pd.notna(row["start"]) and pd.notna(row["end"])
        else pd.NA
    ),
    axis=1
)

# Compute full protein length
df_merged["sequence_length"] = df_merged["Sequence"].str.len()

# Compute domain length
df_merged["domain_length"] = df_merged["domain_sequence"].str.len()
df_merged["domain_length"] = df_merged["domain_length"].astype("Int64")


##### SAVE MERGED TABLE #####
# Save merged table with domain sequences
domain_table_path = os.path.join(results_dir, "0_domain_seq.tsv")
df_merged.to_csv(domain_table_path, sep="\t", index=False)
print("File 0 created at:", domain_table_path)


##### WRITE ONE FASTA PER DOMAIN #####
# Group rows by domain_id and write FASTA files
for domain_id, df_group in df_merged.groupby("domain_id"):

    # Keep only valid domain sequences (exclude NA values)
    df_group = df_group[df_group["domain_sequence"].notna()]
    if df_group.empty:
        continue  # Skip domains without sequences

    # Define FASTA filename and path
    filename = f"{domain_id}.fasta"
    filepath = os.path.join(domains_dir, filename)

    # Write FASTA file
    with open(filepath, "w") as f:
        for _, row in df_group.iterrows():

            uniprot = row["uniprot"]
            start = row["start"]
            end = row["end"]
            seq = row["domain_sequence"]

            # Write FASTA header with UniProt ID, domain ID, and coordinates
            f.write(f">{uniprot}|{domain_id}|{start}-{end}\n")

            # Write domain sequence
            f.write(seq + "\n")

print("FASTAs by domain created successfully")



##### COUNT SEQUENCES PER FASTA #####

# List all FASTA files in output folder
fasta_files = [fn for fn in os.listdir(domains_dir) if fn.endswith(".fasta")]

seq_counts = {}  # Dictionary: {fasta_file: number_of_sequences}

# Count sequences in each FASTA file (count FASTA headers)
for fn in fasta_files:
    filepath = os.path.join(domains_dir, fn)

    # Count FASTA headers (each ">" represents one sequence)
    with open(filepath) as f:
        count = sum(1 for line in f if line.startswith(">"))

    seq_counts[fn] = count

# Convert to DataFrame
df_seq_counts = pd.DataFrame(
    list(seq_counts.items()),
    columns=["fasta_file", "n_sequences"]
)

print(df_seq_counts.head())



##### SAVE SEQUENCE COUNTS #####
counts_path = os.path.join(results_dir, "1_fasta_seq_counts.tsv")
df_seq_counts.to_csv(counts_path, sep="\t", index=False)
print("File 1 created at:", counts_path)



##### SIMPLE DISTRIBUTION (0,1,2,...) #####

# Count how many FASTA files have 0, 1, 2, ... sequences
distribution = df_seq_counts["n_sequences"].value_counts().sort_index()

# Convert Series to DataFrame
df_distribution = distribution.reset_index()
df_distribution.columns = ["n_sequences", "n_fasta_files"]

# Save distribution table
distribution_path = os.path.join(results_dir, "2_fasta_sequence_distribution.tsv")
df_distribution.to_csv(distribution_path, sep="\t", index=False)

print("File 2 created at:", distribution_path)



##### RANGE-BASED DISTRIBUTION #####

# Define numeric bins for sequence ranges
bins = [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 50, 100, 200, 500, 1000, float("inf")]

# Define labels for each range
labels = [
    "0",
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
    "11-50",
    "51-100",
    "101-200",
    "201-500",
    "501-1000",
    "1001+"
]

# Assign each FASTA file to a sequence range
df_seq_counts["seq_range"] = pd.cut(
    df_seq_counts["n_sequences"],
    bins=bins,
    labels=labels
)

# Count number of FASTA files per range
df_range_distribution = (
    df_seq_counts["seq_range"]
    .value_counts()
    .reindex(labels, fill_value=0)
    .reset_index()
)

df_range_distribution.columns = ["seq_range", "n_fasta_files"]

# Save range distribution table
range_path = os.path.join(results_dir, "3_fasta_sequence_distribution_ranges.tsv")
df_range_distribution.to_csv(range_path, sep="\t", index=False)

print("Range distribution file created at:", range_path)