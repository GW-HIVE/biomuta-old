import pandas as pd

# Define file paths
input_file = '/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/mapping_ids/chr_pos_to_ensp.csv'
output_file = '/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/mapping_ids/unique_chr_pos_to_ensp.csv'

# Load the CSV file
df = pd.read_csv(input_file)

# Remove duplicates based on the columns "chr_id", "start_pos", "end_pos", and "entrez"
unique_df = df.drop_duplicates(subset=["chr_id", "start_pos", "end_pos", "entrez"])

# Write the unique rows to a new CSV file
unique_df.to_csv(output_file, index=False)

print(f"Unique entries written to: {output_file}")
