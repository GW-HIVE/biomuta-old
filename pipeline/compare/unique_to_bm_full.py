import pandas as pd

# Load the big and small CSV files
big_file_path = '/data/shared/biomuta/generated/stats/bm_egfr.csv'  # Path to the big CSV file
small_file_path = '/data/shared/biomuta/generated/stats/unique_to_bm.csv'  # Path to the small CSV file
unique_to_bm_full = '/data/shared/biomuta/generated/stats/unique_to_bm_full.csv'

# Read the files into pandas DataFrames
big_df = pd.read_csv(big_file_path)
small_df = pd.read_csv(small_file_path)

# Extract the 7th to 9th columns (columns index 6, 7, and 8 in 0-indexing)
big_subset = big_df.iloc[:, 6:9]

# Ensure small_df has the same columns (7th to 9th columns) for comparison
small_subset = small_df.iloc[:, 0:3]  # small_df column indices

# Filter the rows in the big file where the values in the 7th to 9th columns match any row in the small file by comparing entire rows at once
filtered_df = big_df[big_subset.apply(tuple, axis=1).isin(small_subset.apply(tuple, axis=1))]

# Sort the dataframe by the 7th column (AA position)
sorted_df = filtered_df.sort_values(by=filtered_df.columns[6])

# Write the filtered rows to a new CSV file
sorted_df.to_csv(unique_to_bm_full, index=False)

# Define the tuple of column indices that should be unique (aa_pos, ref_aa, alt_aa)
unique_column_indices = [6, 7, 8]

# Drop duplicate rows based on the unique tuple of column indices
df_unique = sorted_df.drop_duplicates(subset=sorted_df.columns[unique_column_indices])

# Now count instances of unique values in the 12th column (source)
value_counts = df_unique.iloc[:, 11].value_counts()

# Display the result
print(value_counts)

print(f"Filtered and sorted rows have been saved to " + unique_to_bm_full)