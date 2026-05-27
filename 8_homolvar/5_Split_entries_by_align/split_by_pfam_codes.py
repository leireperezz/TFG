# Divide the merged variants file into several files based on the Pfam alignment in which are found.

import os
import pandas as pd
import argparse



# Ask the user which version of gnomAD files is desired
version_files = input("Write the desired gnomAD files version. You can choose between:\n\tfreq-2\n\tfreq-3\n\tfreq-4\n\tfreq-5\n\tfreq-6\n\tfreq-7\n\tfreq-8\n\tno_freq\nPlease, write the version as indicated: ")
version_pattern = '' if version_files == 'no_freq' else '_' + version_files

#full_vs_seed = input("Write the desired alignment type. (Mix condition prioritizes seed alignment and uses full ones when seed is not found)\nYou can choose between:\n\tfull\n\tseed\n\tseed+full(mix)\nPlease, write the version as indicated: ")


# Ask the user which residues computation file type is desired
#case_files = input("Write the file type that you want to use. You can choose between:\n\tupper_case(up)\n\tall_case(all)\nPlease, write the version as indicated between parenthesis: ")
#up_low_case = 'upperCase' if case_files == 'up' else 'keep_pfam'

# Set data path
data_path = "./data"



all_db_entries = f"{data_path}/merged_db_interpro_REC{version_pattern}.txt"

all_db_entries = os.path.normpath(all_db_entries)


# Read merged_db file
merged_db_df = pd.read_csv(all_db_entries, sep="\t", header=0)

# Order the columns in a more readable way
column_order = ['Pfam_found', 'initial_aa', 'final_aa', 'Pos_align', 'position', 'Initial_aa_coincidence', 'database', 'gene', 'NM', 'rsid', 'enst', 'Uniprot_entry', 'Uniprot_entry_name', 'prot_change', 'clinical_significance', 'Condition(s)', 'Review status', 'Accession', 'VariationID', 'AlleleID(s)', 'Canonical SPDI', 'genome_af', 'exome_af']
merged_db_df = merged_db_df[column_order]

# Order by Pfam code
merged_db_df = merged_db_df.sort_values('Pfam_found')

# Group by Pfam code
grouped_db_df = merged_db_df.groupby('Pfam_found')

# Iterate over all the Pfam groups
for pfam, group in grouped_db_df:
    # Sort group by the desired columns
    sorted_group = group.sort_values(['initial_aa', 'final_aa', 'Pos_align', 'Initial_aa_coincidence'])
    
    
    subfolder_path = f"{data_path}/REC_align/pfam_groups_interpro{version_pattern}"
    subfolder_path = os.path.normpath(subfolder_path)

    if not os.path.exists(subfolder_path):
        # If does not exist, create it
        os.makedirs(subfolder_path)
    
    # Set file name
    file_name = f"{subfolder_path}/{pfam}.txt"
    file_name = os.path.normpath(file_name)
    # Save the group to a CSV file
    sorted_group.to_csv(file_name, sep="\t", index=False)
