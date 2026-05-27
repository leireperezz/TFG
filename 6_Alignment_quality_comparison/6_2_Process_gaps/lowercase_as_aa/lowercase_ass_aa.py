# Remove all-gap columns from filtered alignments, treating lowercase letters as amino acids

import os
import pandas as pd

##### DEFINE INPUT/OUTPUT PATHS #####

datasets = {
    "REC_and_SEED_rec": ("../6_1_Filter_shared_sequences/results/rec_and_seed/REC_filtered", "./results/rec_and_seed/REC_filtered"),
    "REC_and_SEED_seed": ("../6_1_Filter_shared_sequences/results/rec_and_seed/SEED_filtered", "./results/rec_and_seed/SEED_filtered"),
    "REC_and_FULL_rec": ("../6_1_Filter_shared_sequences/results/rec_and_full/REC_filtered", "./results/rec_and_full/REC_filtered"),
    "REC_and_FULL_full": ("../6_1_Filter_shared_sequences/results/rec_and_full/FULL_filtered", "./results/rec_and_full/FULL_filtered")
}

##### PROCESS ALL DATASETS #####

for dataset_name, (input_dir, output_dir) in datasets.items():
    
    os.makedirs(output_dir, exist_ok=True)
    
    for filename in os.listdir(input_dir):
        input_path = os.path.join(input_dir, filename)
        output_path = os.path.join(output_dir, filename)
        
        # Read alignment
        df_pos = (
            pd.read_csv(input_path, header=None)
            .iloc[:, 0]
            .apply(list)
            .apply(pd.Series)
        )
        
        # Convert "." to "-"
        df_pos = df_pos.apply(lambda col: col.map(lambda x: "-" if str(x) == "." else str(x)))
        
        # Remove columns that are all gaps
        df_pos = df_pos.loc[:, ~(df_pos == "-").all(axis=0)]
        
        # Reconstruct sequences
        filtered_seqs = df_pos.fillna("").astype(str).agg("".join, axis=1)
        
        # Save to output directory
        filtered_seqs.to_csv(output_path, index=False, header=False)
    
    print(f"Processed {dataset_name}: {len(os.listdir(input_dir))} files")

print("\nAll alignments processed successfully")