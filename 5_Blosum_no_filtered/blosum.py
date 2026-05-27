# Calculate BLOSUM62 scores for REC, SEED, and FULL alignments without filtering to evaluate individual alignment quality

import os
import math
import pandas as pd
import blosum as bl
import itertools

##### PATHS #####

# Define alignment types and their input directories
alignment_types = {
    "REC": "../4_Generate_alignments/4_3_Postprocess_alignments/results",
    "SEED": "./data/Pfam_SEED_HUMAN_2024_03_20",
    "FULL": "./data/Pfam_FULL_HUMAN_2024_03_20"
}

# Output directory for BLOSUM scores
output_dir = "./results"
os.makedirs(output_dir, exist_ok=True)

##### BLOSUM62 MATRIX #####

blosum62 = bl.BLOSUM(62)

##### BLOSUM62 CALCULATION FUNCTION #####

def calc_blosum62(col):
    """
    Calculate BLOSUM62 score for a single alignment column.
    
    Args:
        col (Series): Column of amino acids from alignment
    
    Returns:
        score_sum (int): Total BLOSUM62 score for all pairwise comparisons
        n_comparisons (int): Number of valid pairwise comparisons
        u_detected (bool): Whether selenocysteine (U) was detected in this column
    """
    # Convert to string and clean None/NaN values
    col_list = [str(x) if pd.notna(x) else "-" for x in col]
    
    score_sum = 0
    n_comparisons = 0
    u_detected = False

    # Calculate pairwise BLOSUM62 scores
    for aa1, aa2 in itertools.combinations(col_list, 2):

        # Skip comparisons involving gaps
        if aa1 == "-" or aa2 == "-":
            continue

        # Skip comparisons involving selenocysteine (U)
        if aa1 == "U" or aa2 == "U":
            u_detected = True
            continue

        # Calculate BLOSUM62 score
        try:
            score = blosum62[aa1][aa2]
        except Exception as e:
            print(f"Unexpected error with {aa1}-{aa2}: {e}")
            raise

        score_sum += score
        n_comparisons += 1

    return score_sum, n_comparisons, u_detected

##### PROCESS EACH ALIGNMENT TYPE #####

for align_type, input_path in alignment_types.items():
    
    print(f"\n{'='*50}")
    print(f"Processing {align_type} alignments")
    print(f"{'='*50}")
    
    output_file = os.path.join(output_dir, f"{align_type.lower()}_blosum.csv")
    
    all_results = []

    # Counters for final summary
    pfams_gt1_no_valid_comparisons = []
    pfams_gt1_no_blosum_due_to_u = []

    for filename in os.listdir(input_path):
        file_path = os.path.join(input_path, filename)

        if not os.path.isfile(file_path):
            continue

        # Read alignment file
        try:
            df_pos = (
                pd.read_csv(file_path, header=None, sep=r"\s+", usecols=[1])
                .iloc[:, 0]
                .apply(list)
                .apply(pd.Series)
            )

        except pd.errors.EmptyDataError:
            print(f"Empty file: {filename}")
            continue
        except ValueError:
            print(f"File without valid columns: {filename}")
            continue
        except Exception as e:
            print(f"Skipping: {filename} | {str(e)}")
            continue

        ##### PREPROCESSING #####

        # Fill NaN with gaps
        df_pos = df_pos.fillna("-")
        
        # Convert lowercase and dots to uppercase or gaps
        df_pos = df_pos.apply(
            lambda col: col.map(lambda x: "-" if str(x) == "." else str(x).upper())
        )

        # Remove columns that are all gaps
        df_pos = df_pos.loc[:, ~(df_pos == "-").all(axis=0)]

        # Number of sequences in alignment
        n_sequences_alignment = df_pos.shape[0]

        # Extract Pfam code from filename
        pfam_code = os.path.splitext(filename)[0].split("_")[0]

        ##### CALCULATE BLOSUM62 SCORE #####

        total_score_sum = 0
        total_n_comparisons = 0
        u_problem_in_file = False

        # Calculate BLOSUM62 for each column
        for col_name in df_pos.columns:
            col_score_sum, col_n_comparisons, col_u_problem = calc_blosum62(df_pos[col_name])
            total_score_sum += col_score_sum
            total_n_comparisons += col_n_comparisons
            if col_u_problem:
                u_problem_in_file = True

        # Calculate global BLOSUM62 mean
        if total_n_comparisons == 0:
            global_blosum62 = float("nan")

            if n_sequences_alignment > 1:
                pfams_gt1_no_valid_comparisons.append(pfam_code)
                print(f"{pfam_code} -> >1 seq BUT 0 valid comparisons")

        else:
            global_blosum62 = total_score_sum / total_n_comparisons

            # Check if inf (can happen if only U residues were present)
            if math.isinf(global_blosum62):
                if u_problem_in_file:
                    pfams_gt1_no_blosum_due_to_u.append(pfam_code)
                    print(f"{pfam_code} -> Problem with U (infinite score)")

            global_blosum62 = round(global_blosum62, 2)

        ##### SAVE RESULT #####

        df_result = pd.DataFrame({
            "Pfam_code": [pfam_code],
            "blosum62_global_mean": [global_blosum62],
        })

        all_results.append(df_result)

    ##### SAVE FINAL CSV FILE #####

    if all_results:
        final_df = pd.concat(all_results, ignore_index=True)
        final_df.to_csv(output_file, index=False)
        print(f"\nFile saved at: {output_file}")
    else:
        final_df = pd.DataFrame(columns=["Pfam_code", "blosum62_global_mean"])
        print("\nNo valid files were processed")

    ##### FINAL SUMMARY #####

    # Remove duplicates for safety
    n_gt1_no_valid_comparisons = len(set(pfams_gt1_no_valid_comparisons))
    n_gt1_no_blosum_due_to_u = len(set(pfams_gt1_no_blosum_due_to_u))

    # Count Pfams with valid BLOSUM
    n_valid_blosum = final_df["blosum62_global_mean"].notna().sum()
    n_invalid_blosum = final_df["blosum62_global_mean"].isna().sum()

    print("\n=========================")
    print(f"{align_type} SUMMARY")
    print("=========================")

    print(f"Total Pfams processed: {len(final_df)}")
    print(f"Pfams with valid BLOSUM: {n_valid_blosum}")
    print(f"Pfams without valid BLOSUM (NaN): {n_invalid_blosum}")

    print(f"\nPfams with >1 seq but 0 valid comparisons: {n_gt1_no_valid_comparisons}")
    print(f"Pfams with BLOSUM problem due to U: {n_gt1_no_blosum_due_to_u}")

print("\n" + "="*50)
print("ALL ALIGNMENT TYPES PROCESSED SUCCESSFULLY")
print("="*50)