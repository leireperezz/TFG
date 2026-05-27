# Plot conservation percentage improvement distribution for REC vs FULL alignments

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

##### INPUT/OUTPUT PATHS #####

# Input file containing REC vs FULL conservation comparison
input_path = "../B_Merge_and_comparison/results_merge/conservation_rec_vs_full_comparison.csv"

# Output path for the plot
output_path = "./results_plot/conservation_plot.png"

##### READ INPUT FILE #####

df = pd.read_csv(input_path)

##### BUILD LONG FORMAT TABLE #####

# Keep only REC vs FULL conservation percentage improvement
df_long = pd.DataFrame({
    "comparison": ["REC vs FULL"] * len(df),
    "improvement_percent": df["REC_vs_FULL_conservation_%_improvement"]
})

##### CLEAN DATA #####

# Convert values to numeric and remove NaN values
df_long["improvement_percent"] = pd.to_numeric(df_long["improvement_percent"], errors="coerce")
df_long = df_long.dropna()

##### DEFINE IMPROVEMENT INTERVALS #####

# Define bins for percentage improvement
bins = [-1000000, -100, -50, -40, -30, -20, -10, 0, 0.000001, 10, 20, 30, 40, 50, 100, 1000000]

labels = [
    "< -100",
    "-100 to -50",
    "-50 to -40",
    "-40 to -30",
    "-30 to -20",
    "-20 to -10",
    "-10 to 0",
    "0",
    "0 to 10",
    "10 to 20",
    "20 to 30",
    "30 to 40",
    "40 to 50",
    "50 to 100",
    ">100"
]

# Assign each alignment to an improvement interval
df_long["interval"] = pd.cut(
    df_long["improvement_percent"],
    bins=bins,
    labels=labels,
    right=False,
    include_lowest=True
)

##### COUNT ALIGNMENTS PER INTERVAL #####

# Count how many alignments fall into each interval
df_counts = (
    df_long.groupby("interval", observed=False)
    .size()
    .reset_index(name="count")
)

##### CREATE PLOT #####

sns.set_style("whitegrid")
plt.figure(figsize=(14, 6))

ax = sns.barplot(
    data=df_counts,
    x="interval",
    y="count"
)

# Add value labels on top of each bar
for i, row in df_counts.iterrows():
    plt.text(
        i,
        row["count"] + 1,
        int(row["count"]),
        ha='center'
    )

plt.xlabel("% improvement interval")
plt.ylabel("Number of alignments")
plt.title("Distribution of conservation percentage improvement (REC vs FULL) (lowercase as aa)")
plt.xticks(rotation=45)
plt.tight_layout()

plt.savefig(output_path, dpi=300)
plt.show()

print(f"Conservation plot saved at: {output_path}")