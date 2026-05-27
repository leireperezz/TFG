# Merge and compare quality metrics (conservation, BLOSUM, gaps) between REC and FULL alignments (lowercase as gaps)

import pandas as pd

##### INPUT PATHS #####

# REC and FULL metric files from previous step
rec_path = "../A_Score_and_conservation/results/rec_cons_score_%_gaps_lowercase_as_gap.csv"
full_path = "../A_Score_and_conservation/results/full_cons_score_%_gaps_lowercase_as_gap.csv"

##### OUTPUT PATHS #####

merge_output_path = "./results_merge/merge_rec_vs_full_conservation_and_score.csv"
cons_output_path = "./results_merge/conservation_rec_vs_full_comparison.csv"
blosum_output_path = "./results_merge/blosum_rec_vs_full_comparison.csv"
gaps_output_path = "./results_merge/gaps_rec_vs_full_comparison.csv"

##### READ INPUT FILES #####

df_rec = pd.read_csv(rec_path)
df_full = pd.read_csv(full_path)

##### SELECT AND RENAME COLUMNS #####

df_rec = df_rec[["Pfam_code", "mean_%conservation", "blosum62_global_mean", "perc_gaps"]].rename(
    columns={
        "mean_%conservation": "REC_conservation",
        "blosum62_global_mean": "REC_blosum",
        "perc_gaps": "REC_gaps"
    }
)

df_full = df_full[["Pfam_code", "mean_%conservation", "blosum62_global_mean", "perc_gaps"]].rename(
    columns={
        "mean_%conservation": "FULL_conservation",
        "blosum62_global_mean": "FULL_blosum",
        "perc_gaps": "FULL_gaps"
    }
)

##### MERGE DATAFRAMES #####

df_merged = df_rec.merge(df_full, on="Pfam_code", how="outer")

##### ENSURE NUMERIC TYPES #####

numeric_cols = [
    "REC_conservation", "FULL_conservation",
    "REC_blosum", "FULL_blosum",
    "REC_gaps", "FULL_gaps"
]

df_merged[numeric_cols] = df_merged[numeric_cols].apply(pd.to_numeric, errors="coerce")
df_merged[numeric_cols] = df_merged[numeric_cols].replace([float("inf"), -float("inf")], pd.NA)

##### SAVE MERGED FILE #####

df_merged.to_csv(merge_output_path, index=False, float_format="%.2f")
print(f"Merged file saved at: {merge_output_path}")

##### CONSERVATION COMPARISON TABLE #####

df_cons = pd.DataFrame()

df_cons["Pfam_code"] = df_merged["Pfam_code"]
df_cons["REC_conservation"] = df_merged["REC_conservation"]
df_cons["FULL_conservation"] = df_merged["FULL_conservation"]

df_cons["REC_vs_FULL_conservation_diff"] = (
    df_merged["REC_conservation"] - df_merged["FULL_conservation"]
)

df_cons["REC_vs_FULL_conservation_%_improvement"] = (
    (
        df_merged["REC_conservation"] - df_merged["FULL_conservation"]
    ) / df_merged["FULL_conservation"].where(df_merged["FULL_conservation"] != 0)
) * 100

df_cons = df_cons.round(2)

print("\nConservation comparison table:")
print(df_cons.head())

df_cons.to_csv(cons_output_path, index=False, float_format="%.2f")
print(f"Conservation comparison saved at: {cons_output_path}")

##### BLOSUM COMPARISON TABLE #####

df_blosum = pd.DataFrame()

df_blosum["Pfam_code"] = df_merged["Pfam_code"]
df_blosum["REC_blosum"] = df_merged["REC_blosum"]
df_blosum["FULL_blosum"] = df_merged["FULL_blosum"]

df_blosum["REC_vs_FULL_blosum_diff"] = (
    df_merged["REC_blosum"] - df_merged["FULL_blosum"]
)

full_abs = df_merged["FULL_blosum"].abs()

df_blosum["REC_vs_FULL_blosum_%_improvement"] = (
    (
        df_merged["REC_blosum"] - df_merged["FULL_blosum"]
    ) / full_abs.where(full_abs != 0)
) * 100

numeric_cols = df_blosum.columns.drop("Pfam_code")
df_blosum[numeric_cols] = df_blosum[numeric_cols].apply(pd.to_numeric, errors="coerce")

print("\nBLOSUM comparison table:")
print(df_blosum.head())

df_blosum.to_csv(blosum_output_path, index=False, float_format="%.2f")
print(f"BLOSUM comparison saved at: {blosum_output_path}")

##### GAPS COMPARISON TABLE #####

df_gaps = pd.DataFrame()

df_gaps["Pfam_code"] = df_merged["Pfam_code"]
df_gaps["REC_gaps"] = df_merged["REC_gaps"]
df_gaps["FULL_gaps"] = df_merged["FULL_gaps"]

df_gaps["REC_vs_FULL_gaps_diff"] = (
    df_merged["FULL_gaps"] - df_merged["REC_gaps"]
)

df_gaps["REC_vs_FULL_gaps_%_improvement"] = (
    (
        df_merged["FULL_gaps"] - df_merged["REC_gaps"]
    ) / df_merged["FULL_gaps"].where(df_merged["FULL_gaps"] != 0)
) * 100

numeric_cols = df_gaps.columns.drop("Pfam_code")
df_gaps[numeric_cols] = df_gaps[numeric_cols].apply(pd.to_numeric, errors="coerce")

print("\nGaps comparison table:")
print(df_gaps.head())

df_gaps.to_csv(gaps_output_path, index=False, float_format="%.2f")
print(f"Gaps comparison saved at: {gaps_output_path}")

##### COUNT BETTER ALIGNMENTS #####

# CONSERVATION
rec_better_than_full_cons = (df_cons["REC_vs_FULL_conservation_diff"] > 0).sum()
full_better_than_rec_cons = (df_cons["REC_vs_FULL_conservation_diff"] < 0).sum()
equal_cons = (df_cons["REC_vs_FULL_conservation_diff"] == 0).sum()

print("\n" + "="*50)
print("CONSERVATION COMPARISON")
print("="*50)
print(f"REC better than FULL: {rec_better_than_full_cons}")
print(f"FULL better than REC: {full_better_than_rec_cons}")
print(f"Equal: {equal_cons}")
print(f"% better in conservation: {round((rec_better_than_full_cons / (full_better_than_rec_cons + equal_cons + rec_better_than_full_cons)) * 100, 2)}%")
print(f"Mean % improvement: {round(df_cons['REC_vs_FULL_conservation_%_improvement'].mean(), 2)}%")

# BLOSUM
rec_better_than_full_blosum = (df_blosum["REC_vs_FULL_blosum_diff"] > 0).sum()
full_better_than_rec_blosum = (df_blosum["REC_vs_FULL_blosum_diff"] < 0).sum()
equal_blosum = (df_blosum["REC_vs_FULL_blosum_diff"] == 0).sum()

print("\n" + "="*50)
print("BLOSUM COMPARISON")
print("="*50)
print(f"REC better than FULL: {rec_better_than_full_blosum}")
print(f"FULL better than REC: {full_better_than_rec_blosum}")
print(f"Equal: {equal_blosum}")
print(f"% better in BLOSUM: {round((rec_better_than_full_blosum / (full_better_than_rec_blosum + equal_blosum + rec_better_than_full_blosum)) * 100, 2)}%")
print(f"Mean BLOSUM % improvement: {round(df_blosum['REC_vs_FULL_blosum_%_improvement'].mean(), 2)}%")

# GAPS (lower is better)
rec_better_than_full_gaps = (df_gaps["REC_vs_FULL_gaps_diff"] > 0).sum()
full_better_than_rec_gaps = (df_gaps["REC_vs_FULL_gaps_diff"] < 0).sum()
equal_gaps = (df_gaps["REC_vs_FULL_gaps_diff"] == 0).sum()

print("\n" + "="*50)
print("GAPS COMPARISON (lower is better)")
print("="*50)
print(f"REC better than FULL (fewer gaps): {rec_better_than_full_gaps}")
print(f"FULL better than REC (fewer gaps): {full_better_than_rec_gaps}")
print(f"Equal: {equal_gaps}")
print(f"% better in gaps: {round((rec_better_than_full_gaps / (full_better_than_rec_gaps + equal_gaps + rec_better_than_full_gaps)) * 100, 2)}%")
print(f"Mean gaps % improvement: {round(df_gaps['REC_vs_FULL_gaps_%_improvement'].mean(), 2)}%")

##### TOP 10 WORST PERFORMING PFAMS (REC vs FULL) #####

# CONSERVATION
top10_worst_cons = (
    df_cons.sort_values("REC_vs_FULL_conservation_diff", ascending=True)
           .head(10)
)

print("\n" + "="*50)
print("Top 10 Pfam with worst performance in REC vs FULL (CONSERVATION)")
print("="*50)
print(
    top10_worst_cons[
        ["Pfam_code", "REC_conservation", "FULL_conservation", "REC_vs_FULL_conservation_diff"]
    ].to_string(index=False)
)

# BLOSUM
top10_worst_blosum = (
    df_blosum.sort_values("REC_vs_FULL_blosum_diff", ascending=True)
             .head(10)
)

print("\n" + "="*50)
print("Top 10 Pfam with worst performance in REC vs FULL (BLOSUM)")
print("="*50)
print(
    top10_worst_blosum[
        ["Pfam_code", "REC_blosum", "FULL_blosum", "REC_vs_FULL_blosum_diff"]
    ].to_string(index=False)
)

# GAPS (negative diff = REC has more gaps than FULL = worse)
top10_worst_gaps = (
    df_gaps.sort_values("REC_vs_FULL_gaps_diff", ascending=True)
           .head(10)
)

print("\n" + "="*50)
print("Top 10 Pfam with worst performance in REC vs FULL (GAPS - more gaps)")
print("="*50)
print(
    top10_worst_gaps[
        ["Pfam_code", "REC_gaps", "FULL_gaps", "REC_vs_FULL_gaps_diff"]
    ].to_string(index=False)
)