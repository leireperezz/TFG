# Calculate conservation percentage, BLOSUM62 score, and gap percentage for REC and SEED alignments (lowercase as amino acids)

import os
import pandas as pd
import blosum as bl
import itertools

##### BLOSUM62 MATRIX #####

blosum62 = bl.BLOSUM(62)

##### CONSERVATION CALCULATION FUNCTION #####

def calc_conservation(col):
    """
    Calculate column conservation percentage.
    
    Args:
        col (Series): Column of amino acids from alignment
    
    Returns:
        conservation (float): Percentage of most frequent uppercase amino acid
    """
    col = col.astype(str)

    total_positions = len(col)
    
    # If only 1 sequence, conservation is not meaningful
    if total_positions <= 1:
        return float('nan')
    
    # Filter only uppercase amino acids (gaps and lowercase excluded)
    aa_upper = col[col.str.fullmatch(r"[A-Z]")]

    if len(aa_upper) == 0:
        return 0

    # Calculate conservation as frequency of most common amino acid
    max_count = aa_upper.value_counts().max()
    return max_count / total_positions

##### BLOSUM62 CALCULATION FUNCTION #####

def calc_blosum62(col):
    """
    Calculate BLOSUM62 score for a single alignment column.
    
    Args:
        col (Series): Column of amino acids from alignment
    
    Returns:
        score_sum (int): Total BLOSUM62 score for all pairwise comparisons
        n_comparisons (int): Number of valid pairwise comparisons
    """
    col_list = col.astype(str).tolist()
    score_sum = 0
    n_comparisons = 0

    # Calculate pairwise BLOSUM62 scores
    for aa1, aa2 in itertools.combinations(col_list, 2):

        # Skip comparisons involving gaps
        if aa1 == "-" or aa2 == "-":
            continue

        score = blosum62[aa1][aa2]
        score_sum += score
        n_comparisons += 1

    return score_sum, n_comparisons

##### GAP PERCENTAGE CALCULATION FUNCTION #####

def calc_gap_percentage(df_pos):
    """
    Calculate global gap percentage of the alignment.
    
    Args:
        df_pos (DataFrame): Alignment data (rows = sequences, columns = positions)
    
    Returns:
        gap_percentage (float): Percentage of gaps in the entire alignment
    """
    total_positions = df_pos.shape[0] * df_pos.shape[1]
    gap_count = (df_pos == "-").sum().sum()

    if total_positions == 0:
        return 0

    return round(gap_count / total_positions, 2)

##### DEFINE DATASETS TO PROCESS #####

datasets = {
    "REC": ("../../6_2_Process_gaps/lowercase_as_aa/results/rec_and_seed/REC_filtered", "./results/rec_cons_score_%_gaps_lowercase_as_aa.csv"),
    "SEED": ("../../6_2_Process_gaps/lowercase_as_aa/results/rec_and_seed/SEED_filtered", "./results/seed_cons_score_%_gaps_lowercase_as_aa.csv")
}

##### PROCESS BOTH DATASETS #####

for align_type, (input_path, output_file) in datasets.items():
    
    print(f"\nProcessing {align_type} alignments...")
    
    output_dir = os.path.dirname(output_file)
    os.makedirs(output_dir, exist_ok=True)
    
    all_results = []

    for filename in os.listdir(input_path):
        file_path = os.path.join(input_path, filename)

        try:
            # Read alignment file (each row is one aligned sequence)
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
        
        # Convert lowercase amino acids to uppercase
        df_pos = df_pos.apply(lambda col: col.map(lambda x: str(x).upper()))

        ##### CALCULATE CONSERVATION #####

        conservation = df_pos.apply(calc_conservation, axis=0)
        mean_conservation = round(conservation.mean(), 2)

        ##### CALCULATE BLOSUM62 SCORE #####

        total_score_sum = 0
        total_n_comparisons = 0

        for col_name in df_pos.columns:
            col_score_sum, col_n_comparisons = calc_blosum62(df_pos[col_name])
            total_score_sum += col_score_sum
            total_n_comparisons += col_n_comparisons

        if total_n_comparisons == 0:
            global_blosum62 = float("nan")
        else:
            global_blosum62 = round(total_score_sum / total_n_comparisons, 2)

        ##### CALCULATE GAP PERCENTAGE #####

        perc_gaps = calc_gap_percentage(df_pos)

        ##### EXTRACT PFAM CODE #####

        pfam_code = os.path.splitext(filename)[0].split("_")[0]

        ##### SAVE RESULT #####

        df_result = pd.DataFrame({
            "Pfam_code": [pfam_code],
            "mean_%conservation": [mean_conservation],
            "blosum62_global_mean": [global_blosum62],
            "perc_gaps": [perc_gaps]
        })

        all_results.append(df_result)

    ##### SAVE FINAL CSV FILE #####

    if all_results:
        final_df = pd.concat(all_results, ignore_index=True)
        final_df.to_csv(output_file, index=False)
        print(f"{align_type} file saved at: {output_file}")
    else:
        print(f"\nNo valid {align_type} files were processed")

print("\nAll alignments processed successfully")