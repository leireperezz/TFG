# Query InterPro API to retrieve Pfam domain coordinates for all human reviewed proteins from UniProt

import pandas as pd
import requests
import os


def load_uniprot_ids(file_name, file_path):
    """
    Load UniProt IDs from a TSV file.
    
    Args:
        file_name (str): Name of the TSV file
        file_path (str): Path to the folder containing the file
    
    Returns:
        uniprot_ids (list): Cleaned list of UniProt IDs from the 'Entry' column
    """
    # Build the full file path by joining folder path and file name
    full_path = os.path.join(file_path, file_name)

    # Read TSV file as strings to preserve ID format
    df = pd.read_csv(full_path, sep="\t", dtype=str)

    # Ensure required column exists
    if "Entry" not in df.columns:
        raise ValueError("The file does not contain the 'Entry' column.")
    
    # Clean IDs: remove NA values, strip whitespace, and remove duplicates
    uniprot_ids = (
        df["Entry"]
        .dropna()
        .str.strip()
        .drop_duplicates()
        .tolist()
    )

    return uniprot_ids


def fetch_pfam_domains_to_file(uniprot_id, file_handle, page_size=200):
    """
    Query InterPro API for Pfam domains of a given UniProt ID and write results directly to an open TSV file.
    
    Args:
        uniprot_id (str): UniProt accession ID
        file_handle (file object): Open file handle for writing TSV output
        page_size (int): Number of results per API page (default: 200)
    
    Returns:
        None (writes directly to file)
    """   
    # Ensure UniProt ID is a clean string (remove surrounding spaces)
    uniprot_id = str(uniprot_id).strip()

    # Build InterPro API URL for Pfam domain queries
    url = (
        f"https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/UniProt/{uniprot_id}/"
        f"?extra_fields=short_name&page_size={page_size}"
    )

    wrote_any = False      # Tracks whether at least one domain row has been written for this protein
    next_url = url         # Used for API pagination

    try:
        while next_url:    # Continue while there is a URL to request (handles pagination)
            # Send GET request to the InterPro API
            response = requests.get(next_url, headers={"Accept": "application/json"})

            if response.status_code == 200:
                data = response.json()           # Convert JSON response into a Python dictionary
                next_url = data.get("next")      # Get next page URL (if available for pagination)

                # Iterate over Pfam domain entries in the API response
                for entry in data.get("results", []): 
                    domain_id = entry["metadata"]["accession"]
                    domain_name = entry["metadata"]["name"]

                    # Extract domain coordinates for each protein match
                    for protein in entry.get("proteins", []):
                        for loc in protein.get("entry_protein_locations", []):
                            for fragment in loc.get("fragments", []):
                                start = fragment.get("start")
                                end = fragment.get("end")

                                # Skip fragments without valid start/end coordinates
                                if start is None or end is None:
                                    continue    

                                # Write domain information to the TSV file
                                file_handle.write(
                                    f"{uniprot_id}\t{domain_id}\t{domain_name}\t{start}\t{end}\n"
                                )
                                wrote_any = True

            elif response.status_code == 204:
                break         # No content returned (protein has no Pfam domains)
            else:
                break         # HTTP error occurred

    except Exception:
        pass                  # Ignore network or runtime errors and continue with next protein

    # If no domains were written, add a row with NA values to keep protein in the output
    if not wrote_any:
            file_handle.write(f"{uniprot_id}\tNA\tNA\tNA\tNA\n")



##### MAIN SCRIPT #####

# Step 1: Load UniProt IDs from input file
input_folder = "../1_Table_uniprot_data"
input_folder = os.path.normpath(input_folder)
input_filename = "uniprotkb_organism_id_9606_AND_reviewed_2026_02_04.tsv"

# Load UniProt IDs from the TSV file
ids = load_uniprot_ids(input_filename, input_folder)

print(f"Loaded {len(ids)} proteins from UniProt")
print(f"First 10 UniProt IDs: {ids[:10]}")

# Step 2: Query InterPro API and create output TSV file with Pfam domain coordinates
output_file = "./pfam_domains_all.tsv"
output_file = os.path.normpath(output_file)

# Create (overwrite) the output TSV file
with open(output_file, "w", encoding="utf-8") as f:
    
    # Write TSV header
    f.write("uniprot\tdomain_id\tdomain_name\tstart\tend\n")

    # Query the API for each UniProt ID and write rows to the file
    for uniprot_id in ids:
        fetch_pfam_domains_to_file(uniprot_id, f)

# Print confirmation message with output file location
print(f"Output file created: {output_file}")
