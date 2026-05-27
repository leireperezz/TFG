# Plot BLOSUM62 score comparison distributions (percentage improvement and absolute difference) for REC vs SEED alignments (lowercase as gaps)

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

##### INPUT/OUTPUT PATHS #####

# Input file containing REC vs SEED BLOSUM comparison
input_path = "../B_Merge_and_comparison/results_merge/blosum_rec_vs_seed_comparison.csv"

# Output paths for plots
percent_output_path = "./results_plot/blosum_percent_plot.png"
diff_output_path = "./results_plot/blosum_diff_plot.png"

##### READ INPUT FILE #####

df = pd.read_csv(input_path)

##### GENERAL PLOT STYLE #####

sns.set_style("whitegrid")
sns.set_palette("Set2")

##### PLOT 1: BLOSUM % IMPROVEMENT DISTRIBUTION #####

# Prepare data for percentage improvement plot
df_percent_long = pd.DataFrame({
    "comparison": ["REC vs SEED"] * len(df),
    "value": df["REC_vs_SEED_blosum_%_improvement"]
})

# Convert values to numeric and remove missing values
df_percent_long["value"] = pd.to_numeric(df_percent_long["value"], errors="coerce")
df_percent_long = df_percent_long.dropna()

# Define bins for percentage improvement intervals
bins_percent = [
    -1000000, -200, -100, -50, -25, -10, -5,
    0, 0.000001,
    5, 10, 25, 50, 100, 200, 1000000
]
labels_percent = [
    "< -200", "-200 to -100", "-100 to -50", "-50 to -25", "-25 to -10", "-10 to -5",
    "-5 to 0", "0", "0 to 5", "5 to 10", "10 to 25", "25 to 50", "50 to 100",
    "100 to 200", "> 200"
]

# Assign each alignment to a percentage improvement interval
df_percent_long["interval"] = pd.cut(
    df_percent_long["value"],
    bins=bins_percent,
    labels=labels_percent,
    right=False,
    include_lowest=True
)

# Count alignments per interval
df_percent_counts = (
    df_percent_long.groupby("interval", observed=False)
    .size()
    .reset_index(name="count")
)

# Create percentage improvement plot
plt.figure(figsize=(14, 6))

ax = sns.barplot(
    data=df_percent_counts,
    x="interval",
    y="count"
)

# Add value labels on top of bars
for i, row in df_percent_counts.iterrows():
    plt.text(
        i,
        row["count"] + 1,
        int(row["count"]),
        ha="center"
    )

plt.xlabel("% improvement interval")
plt.ylabel("Number of alignments")
plt.title("Distribution of percentage change in BLOSUM62 score (REC vs SEED) (lowercase as gap)")
plt.xticks(rotation=45)
plt.tight_layout()

plt.savefig(percent_output_path, dpi=300)
plt.show()

print(f"BLOSUM % improvement plot saved at: {percent_output_path}")

##### PLOT 2: BLOSUM ABSOLUTE DIFFERENCE DISTRIBUTION #####

# Prepare data for absolute difference plot
df_diff_long = pd.DataFrame({
    "comparison": ["REC vs SEED"] * len(df),
    "value": df["REC_vs_SEED_blosum_diff"]
})

# Convert values to numeric and remove missing or infinite values
df_diff_long["value"] = pd.to_numeric(df_diff_long["value"], errors="coerce")
df_diff_long = df_diff_long.replace([np.inf, -np.inf], np.nan).dropna()

# Define bins for absolute difference intervals
bins_diff = [
    -1000000, -500, -300, -200, -100, -50, -20, -10,
    0, 0.000001,
    10, 20, 50, 100, 200, 300, 500, 1000000
]

labels_diff = [
    "< -500", "-500 to -300", "-300 to -200", "-200 to -100", "-100 to -50", "-50 to -20", "-20 to -10",
    "-10 to 0", "0", "0 to 10", "10 to 20", "20 to 50", "50 to 100", "100 to 200",
    "200 to 300", "300 to 500", "> 500"
]

# Assign each alignment to a difference interval
df_diff_long["interval"] = pd.cut(
    df_diff_long["value"],
    bins=bins_diff,
    labels=labels_diff,
    right=False,
    include_lowest=True
)

# Count alignments per interval
df_diff_counts = (
    df_diff_long.groupby("interval", observed=False)
    .size()
    .reset_index(name="count")
)

# Create absolute difference plot
plt.figure(figsize=(14, 6))

ax = sns.barplot(
    data=df_diff_counts,
    x="interval",
    y="count"
)

# Add value labels on top of bars
for i, row in df_diff_counts.iterrows():
    plt.text(
        i,
        row["count"] + 1,
        int(row["count"]),
        ha="center"
    )

plt.xlabel("BLOSUM62 score difference interval")
plt.ylabel("Number of alignments")
plt.title("Distribution of BLOSUM62 score differences (REC vs SEED) (lowercase as gap)")
plt.xticks(rotation=45)
plt.tight_layout()

plt.savefig(diff_output_path, dpi=300)
plt.show()

print(f"BLOSUM absolute difference plot saved at: {diff_output_path}")