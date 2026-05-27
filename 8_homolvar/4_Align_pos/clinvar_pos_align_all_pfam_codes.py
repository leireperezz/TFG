#tengo que quitar la pregunta de que alineamiento quieres y añadir la ruta de la carpeta con los alineamientos seleccionados
#hay que modificar el nombre del archivo que espera ya que no son iguales

# Get equivalent positions for ClinVar variants using Pfam alignments
# This code was adapted from Sergi Soldevila in order to consider all Pfam codes, not only the first.

import contextlib
import os
import pandas as pd
from Bio.Data.IUPACData import protein_letters_3to1
import glob

def pos_align(clinvar_df):
   
    pfam_path = "../../4_Generate_alignments/4_3_Postprocess_alignments/results"
    
    
    pfam_path = os.path.normpath(pfam_path)

    # Ask the user if wants to compute only upper case residues from Pfam or turn all of them to upper case
    #case_files = input("Pfam alignments have lower case residues indicating less conservation. Write the residues that you want to compute. You can change all residues to upper case (up) or keep Pfam structure (pfam):\n\tupper_case(up)\n\tkeep_pfam(pfam)\nPlease, write the version as indicated between parenthesis: ")
    #up_low_case = 'upperCase' if case_files == 'up' else 'keep_pfam'
    up_low_case = 'keep_pfam'

    # Set a boolean to indicate if the position is already found in the row to avoid checking more Pfam alignments
    switch = False
    # Set a counter to check the cases for which the alignment was not found
    cases = 0

    # Iterate over the rows of ClinVar DataFrame
    for index, row in clinvar_df.iterrows():
        # Save the list of Pfam codes in the current row
        pfam_list = row["Pfam"]
        # Iterate over the list of Pfam codes
        for pfam_code in pfam_list:
            # If the position has not been retrieved for this row yet
            if switch == False:
                # Handle the possibility that the corresponding Pfam alignment file is not found
                with contextlib.suppress(FileNotFoundError):
                    # Set the name of the file from the Pfam code
                    '''
                    pfam_file = pfam_path + "/" + str(pfam_code) + "_HUMAN.txt"
                    pfam_file = os.path.normpath(pfam_file) '''

                    pattern = os.path.join(pfam_path, str(pfam_code) + "_*.txt")
                    files = glob.glob(pattern)

                    if len(files) == 0:
                        continue

                    pfam_file = os.path.normpath(files[0])


                    # print(pfam_file)

                    # Open the corresponding Pfam alignment file
                    with open(pfam_file) as f:
                        # Save the content of the file
                        file_content = f.readlines()
                    # Iterate over the contents of the file
                    for text in file_content:
                        # Remove the newline characters
                        text = text.replace("\n", "")
                        # Save the UniProt entry code and transform it to string
                        uniprot = row["Uniprot_entry_name"]
                        uniprot = str(uniprot)
                        # Check if the text starts by the UniProt code
                        if text.startswith(uniprot):
                            # If yes, save the rang, initial and final position variables from the text
                            rang = text.split("/")[1].split("-", maxsplit=1)
                            pos_i = float(rang[0])
                            pos_f = float(rang[1].split()[0])
                            # Save the number indicating the position of the aminoacidic change and transform it to float
                            aa_change_pos = row["Position"]
                            aa_change_pos = float(aa_change_pos)
                            # Check if the aminoacidic change position fits inside the range of the alignment position
                            if aa_change_pos >= pos_i and aa_change_pos <= pos_f:
                                # If yes, initialize two counters to 0.
                                compt_aa = 0
                                compt = 0
                                # Save the part of the text corresponding to the alignment
                                alignment = text.split()[1]

                                # Iterate over each position of the alignment
                                for align_pos in alignment:
                                    # Add one to the counter each time
                                    compt = compt + 1
                                    # If the position is an aminoacid, add one to the other counter
                                    if align_pos.isalpha() == True:
                                        compt_aa = compt_aa + 1
                                        # Check if the counter of aminoacids is equal to this expression: position of the aminoacid change - initial position in the pfam alignment + 1
                                        if compt_aa == (aa_change_pos - pos_i + 1):
                                            # If yes, save the Pfam code where the alignment was found
                                            clinvar_df.at[index, "Pfam_found"] = pfam_code
                                            # Add the counter of positions in the alignment to the Pos_align column in the corresponding row
                                            clinvar_df.at[index, "Pos_align"] = compt
                                            
                                            # Extract the initial aa
                                            initial_aa = row['Initial aa']
                                            if initial_aa in protein_letters_3to1:
                                                initial_aa = protein_letters_3to1[initial_aa]
                                            # Check if the aminoacid in the alignment is the same as the initial one, compute it depending on the chosen case (considering only uppercase or all residues)
                                            if (up_low_case == 'upperCase' and align_pos.upper() == initial_aa
                                                or up_low_case == 'keep_pfam' and align_pos == initial_aa):
                                                # If yes, introduce True to Initial_aa_coincidence column
                                                clinvar_df.at[index, "Initial_aa_coincidence"] = True
                                            else:
                                                # If not, introduce False
                                                clinvar_df.at[index, "Initial_aa_coincidence"] = False

                                            # change the switch boolean to True, meaning that it was found the alignment and the position is added
                                            switch = True
                                            # Break the for loop
                                            break
                                        # print(compt_aa-1, compt, len(alignment))

                                # Reinitialize both counters
                                compt = 0
                                compt_aa = 0

        # Check if the switch boolean is still false, meaning that after checking all pfam codes, the alignment was not found and any position was added
        if switch == False:
            # If yes, add 0 as the position of alignment
            clinvar_df.at[index, "Pos_align"] = 0
            # Add one to the cases counter, that indicates the amount of cases for which the alignment was not found
            cases += 1
        # Set the switch boolean to False, before start checking the next row
        switch = False


    # print(clinvar_df[["Position", "Pos_align"]])

    print(f"The selected alignment was not found in {cases} cases.")
    
    return clinvar_df



#full_vs_seed = input("Write the desired alignment type. (Mix condition prioritizes seed alignment and uses full ones when seed is not found)\nYou can choose between:\n\tfull\n\tseed\n\tseed+full(mix)\nPlease, write the version as indicated: ")

# Set ClinVar file path
clinvar_path = "../3_Data_filtering/data/clinvar_AD_missense_genes_min3.txt"
clinvar_path = os.path.normpath(clinvar_path)

# Read ClinVar file
clinvar_df = pd.read_csv(clinvar_path, sep="\t", header=0)

# Split the different pfam codes and discart the last position, which is always empty
clinvar_df["Pfam"] = clinvar_df["Pfam"].str.split(";").str[:-1]

# Filter the DataFrame avoiding empty values in Pfam column
clinvar_df = clinvar_df[clinvar_df["Pfam"].notna()]


clinvar_df = pos_align(clinvar_df)

# Remove 0s
filtered_clinvar_df = clinvar_df[clinvar_df["Pos_align"] != 0]

# Check if data directory exists
#output_path = "./data"
output_path = ".\data\REC"

if not os.path.exists(output_path):
    # If does not exist, create it
    os.makedirs(output_path)



output_file = f"{output_path}/clinvar_pos_align_REC.txt"
output_file = os.path.normpath(output_file)

# Save the filtered ClinVar DataFrame to a CSV file
filtered_clinvar_df.to_csv(output_file, header=True, index=False, sep="\t", mode="w")

print(f"Variants with alignment position found: {len(filtered_clinvar_df)}")
print(f"Variants without alignment position: {len(clinvar_df) - len(filtered_clinvar_df)}")
print(f"Total ClinVar variants processed: {len(clinvar_df)}")
