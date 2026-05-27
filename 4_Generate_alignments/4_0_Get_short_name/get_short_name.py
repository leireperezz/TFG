# Query InterPro API to retrieve Pfam short names for each domain code

import pandas as pd
import requests
import os


def get_short_name(pfam_code):
    """
    Retrieve the Pfam short_name from the InterPro API.
    
    Args:
        pfam_code (str): Pfam domain accession code
    
    Returns:
        short_name (str): Short name of the Pfam domain, or None if not available
    """
    url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{pfam_code}?extra_fields=short_name"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            return data.get("extra_fields", {}).get("short_name", None)
        else:
            return None
    except Exception:
        return None



##### LOAD INPUT FILE #####

# Input file containing FASTA filenames from previous step
file_path = "../3_Extract_sequences/results/1_fasta_seq_counts.tsv"
file_path = os.path.normpath(file_path)

# Read TSV file
df = pd.read_csv(file_path, sep="\t")



##### EXTRACT PFAM CODES #####

# Create a new DataFrame with only the fasta_file column
df_new = df[["fasta_file"]].copy()

# Remove file extension by splitting at "." and keeping the first part
df_new["fasta_file"] = df_new["fasta_file"].str.split(".").str[0]

# Rename column to pfam_code
df_new = df_new.rename(columns={"fasta_file": "pfam_code"})

print(df_new.head())



##### PREPARE OUTPUT FILE #####

# Create results folder in the script directory (if it does not exist)
script_dir = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.join(script_dir, "results")
os.makedirs(output_dir, exist_ok=True)

# Define output file path
output_file = os.path.join(output_dir, "pfam_code_short_name.tsv")

# Create output file and write header
with open(output_file, "w") as f:
    f.write("pfam_code\tshort_name\n")



##### QUERY API AND WRITE RESULTS #####

# Iterate over Pfam codes and retrieve short_name from InterPro API
for pfam_code in df_new["pfam_code"]:

    short_name = get_short_name(pfam_code)

    # Append result to file
    with open(output_file, "a") as f:
        f.write(f"{pfam_code}\t{short_name}\n")

print("Process finished.")
print(f"Output file saved at: {output_file}")