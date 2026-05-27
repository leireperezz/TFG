# Select best alignments based on BLOSUM62 score comparisons and copy selected alignment files

import pandas as pd
import os

##### INPUT PATHS - COMPARISON FILES #####

# Comparison files from REC vs SEED and REC vs FULL (lowercase as aa)
input_rec_seed_file = "../6_Alignment_quality_comparison/6_3_Calculate_quality_metrics/lowercase_as_aa/rec_vs_seed/B_Merge_and_comparison/results_merge/blosum_rec_vs_seed_comparison.csv"
input_rec_full_file = "../6_Alignment_quality_comparison/6_3_Calculate_quality_metrics/lowercase_as_aa/rec_vs_full/B_Merge_and_comparison/results_merge/blosum_rec_vs_full_comparison.csv"
##### INPUT PATHS - ORIGINAL UNFILTERED BLOSUM SCORES #####

input_rec_original = "../5_Blosum_no_filtered/results/rec_blosum.csv"
input_seed_original = "../5_Blosum_no_filtered/results/seed_blosum.csv"
input_full_original = "../5_Blosum_no_filtered/results/full_blosum.csv"

##### OUTPUT PATHS #####

output_dir = "./results"
output_file = os.path.join(output_dir, "selected_alignments.csv")
os.makedirs(output_dir, exist_ok=True)

##### ALIGNMENT SOURCE DIRECTORIES #####

# Input directories for each alignment type
seed_dir = "../6_Alignment_quality_comparison/6_1_Filter_shared_sequences/data/Pfam_SEED_HUMAN_2024_03_20"
full_dir = "../6_Alignment_quality_comparison/6_1_Filter_shared_sequences/data/Pfam_HUMAN_InterPro_2023-07-03"
rec_dir = "../4_Generate_alignments/4_3_Postprocess_alignments/results"

output_alignments_dir = os.path.join(output_dir, "selected_alignments")
os.makedirs(output_alignments_dir, exist_ok=True)

##### LOAD COMPARISON DATA #####

df_rec_seed = pd.read_csv(input_rec_seed_file)
df_rec_full = pd.read_csv(input_rec_full_file)

##### LOAD ORIGINAL BLOSUM SCORES #####

df_rec_original = pd.read_csv(input_rec_original)
df_seed_original = pd.read_csv(input_seed_original)
df_full_original = pd.read_csv(input_full_original)

df_rec_original = df_rec_original[["Pfam_code", "blosum62_global_mean"]].rename(
    columns={"blosum62_global_mean": "REC_original_score"}
)

df_seed_original = df_seed_original[["Pfam_code", "blosum62_global_mean"]].rename(
    columns={"blosum62_global_mean": "SEED_original_score"}
)

df_full_original = df_full_original[["Pfam_code", "blosum62_global_mean"]].rename(
    columns={"blosum62_global_mean": "FULL_original_score"}
)

##### PART 1: REC VS SEED SELECTION #####

# Select SEED if better than REC, otherwise keep REC
df_seed_better = df_rec_seed[df_rec_seed["REC_vs_SEED_blosum_diff"] < 0].copy()
df_rec_kept_seed = df_rec_seed[df_rec_seed["REC_vs_SEED_blosum_diff"] >= 0].copy()

output_seed = pd.DataFrame({
    "Pfam_code": df_seed_better["Pfam_code"],
    "alignment_type": "SEED"
})

output_rec_seed = pd.DataFrame({
    "Pfam_code": df_rec_kept_seed["Pfam_code"],
    "alignment_type": "REC"
})

output_df = pd.concat([output_seed, output_rec_seed], ignore_index=True)

##### PART 2: REC VS FULL SELECTION #####

existing_pfams = set(output_df["Pfam_code"])

df_full_repeated = df_rec_full[df_rec_full["Pfam_code"].isin(existing_pfams)].copy()
df_full_new = df_rec_full[~df_rec_full["Pfam_code"].isin(existing_pfams)].copy()

df_full_repeated = df_full_repeated.merge(
    output_df[["Pfam_code", "alignment_type"]],
    on="Pfam_code",
    how="left"
)

df_replace_with_full = df_full_repeated[
    (df_full_repeated["REC_vs_FULL_blosum_diff"] < 0) &
    (df_full_repeated["alignment_type"] == "REC")
].copy()

pfams_to_replace = set(df_replace_with_full["Pfam_code"])

output_df = output_df[~output_df["Pfam_code"].isin(pfams_to_replace)].copy()

output_full_replacements = pd.DataFrame({
    "Pfam_code": df_replace_with_full["Pfam_code"],
    "alignment_type": "FULL"
})

df_full_better_new = df_full_new[df_full_new["REC_vs_FULL_blosum_diff"] < 0].copy()
df_rec_kept_full_new = df_full_new[df_full_new["REC_vs_FULL_blosum_diff"] >= 0].copy()

output_full_new = pd.DataFrame({
    "Pfam_code": df_full_better_new["Pfam_code"],
    "alignment_type": "FULL"
})

output_rec_new = pd.DataFrame({
    "Pfam_code": df_rec_kept_full_new["Pfam_code"],
    "alignment_type": "REC"
})

output_df = pd.concat(
    [output_df, output_full_replacements, output_full_new, output_rec_new],
    ignore_index=True
)

##### PART 3: ADD REC ONLY PFAMS #####

selected_pfams = set(output_df["Pfam_code"])
rec_pfams = set(df_rec_original["Pfam_code"])

missing_rec_pfams = rec_pfams - selected_pfams

df_rec_missing = df_rec_original[
    df_rec_original["Pfam_code"].isin(missing_rec_pfams)
].copy()

df_rec_missing = df_rec_missing[
    df_rec_missing["REC_original_score"].notna()
].copy()

output_rec_missing = pd.DataFrame({
    "Pfam_code": df_rec_missing["Pfam_code"],
    "alignment_type": "REC"
})

output_df = pd.concat([output_df, output_rec_missing], ignore_index=True)

##### FINAL CLEANUP #####

output_df = output_df.sort_values("Pfam_code")
output_df = output_df.drop_duplicates(subset=["Pfam_code"], keep="first")

##### ADD ORIGINAL SCORES #####

output_df = output_df.merge(df_seed_original, on="Pfam_code", how="left")
output_df = output_df.merge(df_rec_original, on="Pfam_code", how="left")
output_df = output_df.merge(df_full_original, on="Pfam_code", how="left")

output_df["score_original"] = output_df.apply(
    lambda row: row["SEED_original_score"] if row["alignment_type"] == "SEED"
    else row["REC_original_score"] if row["alignment_type"] == "REC"
    else row["FULL_original_score"] if row["alignment_type"] == "FULL"
    else pd.NA,
    axis=1
)

output_df = output_df[["Pfam_code", "alignment_type", "score_original"]]

output_df.to_csv(output_file, index=False)

##### FUNCTION TO INDEX ALIGNMENT FILES #####

def index_txt_files(folder):
    """
    Create a dictionary mapping Pfam codes to alignment file paths.
    
    Args:
        folder (str): Path to folder containing alignment files
    
    Returns:
        file_index (dict): Dictionary with Pfam code as key and file path as value
    """
    file_index = {}

    for filename in os.listdir(folder):
        if filename.endswith(".txt") and filename.startswith("PF"):
            pfam_code = filename.split("_")[0].split(".")[0]
            file_index[pfam_code] = os.path.join(folder, filename)

    return file_index

##### INDEX ALIGNMENT FOLDERS #####

seed_files_index = index_txt_files(seed_dir)
full_files_index = index_txt_files(full_dir)
rec_files_index = index_txt_files(rec_dir)

alignment_indexes = {
    "SEED": seed_files_index,
    "FULL": full_files_index,
    "REC": rec_files_index
}

##### COPY SELECTED ALIGNMENTS #####

not_found_list = []

for _, row in output_df.iterrows():

    pfam_code = str(row["Pfam_code"]).strip()
    alignment_type = str(row["alignment_type"]).strip().upper()

    file_index = alignment_indexes[alignment_type]

    if pfam_code not in file_index:
        not_found_list.append(pfam_code)
        continue

    source_file = file_index[pfam_code]

    output_alignment_file = os.path.join(
        output_alignments_dir,
        f"{pfam_code}_{alignment_type}.txt"
    )

    # Copy file using binary mode
    with open(source_file, "rb") as source:
        content = source.read()

    with open(output_alignment_file, "wb") as destination:
        destination.write(content)

##### COUNT REAL FILES IN OUTPUT FOLDER #####

seed_count = 0
full_count = 0
rec_count = 0

for filename in os.listdir(output_alignments_dir):
    if filename.endswith("_SEED.txt"):
        seed_count += 1
    elif filename.endswith("_FULL.txt"):
        full_count += 1
    elif filename.endswith("_REC.txt"):
        rec_count += 1

total_count = seed_count + full_count + rec_count

##### PRINT SUMMARY #####

print("\n" + "="*50)
print("ALIGNMENT SELECTION SUMMARY")
print("="*50)
print(f"Final Pfams in CSV: {len(output_df)}")
print(f"CSV saved at: {output_file}")

print("\n" + "="*50)
print("ALIGNMENT TYPE DISTRIBUTION IN CSV")
print("="*50)
print(f"SEED: {len(output_df[output_df['alignment_type'] == 'SEED'])}")
print(f"FULL: {len(output_df[output_df['alignment_type'] == 'FULL'])}")
print(f"REC: {len(output_df[output_df['alignment_type'] == 'REC'])}")

print("\n" + "="*50)
print("FILES IN OUTPUT FOLDER")
print("="*50)
print(f"Total files: {total_count}")
print(f"SEED files: {seed_count}")
print(f"FULL files: {full_count}")
print(f"REC files: {rec_count}")

if not_found_list:
    print(f"\nAlignments not found ({len(not_found_list)}):")
    for pfam in not_found_list:
        print(f"  {pfam}")

print("\nProcess finished.")