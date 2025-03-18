# Removes rows that contain strings of nucleotides in refNt or altNt instead of having one.
# Removes rows with NAs.
# Substitutes DO IDs.
# Prints all unique values from columns where it makes sense (for example, check the number of unique nucleotides, should be 4; the number of unique aa should be 20 etc.)

import csv
import numpy as np
import pandas as pd

bm_path = "/data/shared/repos/biomuta-old/generated_datasets/compiled/biomuta_v6.csv"

bm = pd.read_csv(bm_path)

# Drop rows with missing data
print(f"Number of rows before dropping NAs: {len(bm.index)}")
bm = bm.replace(["None", "nan", "NaN", "NULL", ""], np.nan)
bm_notna = bm.dropna(axis=0, subset=bm.columns.difference(["sample_name", "dbsnp_id"]))
print(f"Number of rows after dropping NAs: {len(bm_notna.index)}") # Number of rows after dropping NAs
print(f"Unique chromosome IDs are: {bm_notna['chr_id'].unique()}")

# Drop rows with invalid values
## Validate nucleotides
valid_nucleotides = ["A", "T", "G", "C"]
is_valid_nucleotide = bm_notna["ref_nt"].isin(valid_nucleotides) & bm_notna["alt_nt"].isin(valid_nucleotides)
## Validate aminoacids
bm_notna = bm_notna.copy()
bm_notna.loc[:, "ref_aa"] = bm_notna["ref_aa"].str.upper()
bm_notna.loc[:, "alt_aa"] = bm_notna["alt_aa"].str.upper()
mask_digits = bm_notna["ref_aa"].str.isdigit() | bm_notna["alt_aa"].str.isdigit()
mask_invalid_aa = (bm_notna["ref_aa"] == "=") | (bm_notna["alt_aa"].isin(["?", "O", "="]))
is_single_aa = (bm_notna["ref_aa"].str.len() == 1) & (bm_notna["alt_aa"].str.len() == 1)
## Filter
valid_rows = is_valid_nucleotide & is_single_aa & ~mask_digits & ~mask_invalid_aa
bm_filtered = bm_notna[valid_rows]
print(f"Number of rows after dropping strings of nucleotides and invalid AA: {len(bm_filtered.index)}")

# Drop rows with new DO IDs
# List approved DO IDs and delete rows that don't contain these numbers in the "do_name" column
# List of approved DO IDs
approved_doids = {
    "DOID:11054",
    "DOID:1612",
    "DOID:9256",
    "DOID:5041",
    "DOID:11934",
    "DOID:263",
    "DOID:3571",
    "DOID:1324",
    "DOID:10283",
    "DOID:10534",
    "DOID:1781",
    "DOID:363",
    "DOID:4362",
    "DOID:1319",
    "DOID:2531",
    "DOID:3953",
    "DOID:1793",
    "DOID:2394",
    "DOID:4159"
}

# Identify rows with the combined DO ID
combined_mask = bm_filtered["do_name"] == "DOID:5041 / esophageal cancer & DOID:10534 / stomach cancer"
combined_rows = bm_filtered[combined_mask].copy()

# Create two separate sets of rows with different DO names
esophageal_rows = combined_rows.copy()
esophageal_rows["do_name"] = "DOID:5041 / esophageal cancer"
stomach_rows = combined_rows.copy()
stomach_rows["do_name"] = "DOID:10534 / stomach cancer"

# Remove the original combined rows and add the new separate rows
df_result = pd.concat([
    bm_filtered[~combined_mask], # Original data excluding the combined rows
    esophageal_rows,    # New rows for esophageal cancer
    stomach_rows        # New rows for stomach cancer
])

# Replace "DOID:3571 / liver cancer & DOID:4606 / bile duct cancer" with "DOID:3571 / liver cancer" (it turns out, Combined Hepatocellular and Intrahepatic Cholangiocarcinoma is a liver cancer)
df_result = df_result.replace(to_replace={
    "DOID:3571 / liver cancer & DOID:4606 / bile duct cancer": "DOID:3571 / liver cancer",
    "DOID:1781 / thyroid cancer": "DOID:1781 / thyroid gland cancer"
    })
df_result["chr_id"] = df_result["chr_id"].astype(str).str.replace("X", "23")

# Extract DOID part
df_result["extracted_doid"] = df_result["do_name"].str.extract(r"(DOID:\d+)")

# Get the deleted DOIDs
deleted_doids = df_result.loc[~df_result["extracted_doid"].isin(approved_doids), "extracted_doid"].dropna().unique()
'''
# Print deleted DOIDs
print("Deleted DOIDs:", deleted_doids)
'''

# Reset index to have clean indices
df_result = df_result.reset_index(drop=True)

# Keep only approved DOIDs
df_filtered = df_result[df_result["extracted_doid"].isin(approved_doids)].drop(columns=["extracted_doid"])

# Print unique do_name values from the filtered DataFrame
'''
print("\nUnique do_name values (Approved DOIDs only):")
print(df_result.loc[df_result['extracted_doid'].isin(approved_doids), 'do_name'].unique())
'''
print(f"dtype of start_pos: {df_filtered['start_pos'].dtype}")
print(f"dtype of end_pos: {df_filtered['end_pos'].dtype}")
print(f"dtype of aa_pos: {df_filtered['aa_pos'].dtype}")
print(f"Unique chromosome IDs are: {df_filtered['chr_id'].unique()}")
print(f"Unique reference nucleotides are: {df_filtered['ref_nt'].unique()}")
print(f"Unique alt nucleotides are: {df_filtered['alt_nt'].unique()}")
print(f"Unique reference aa are: {df_filtered['ref_aa'].unique()}")
print(f"Unique alt aa are: {df_filtered['alt_aa'].unique()}")
print(f"The dimensionality of the final dataframe is {df_filtered.shape}")

# Save the filtered dataframe to a file
outfile = "/data/shared/repos/biomuta-old/generated_datasets/compiled/biomuta_v6.1.csv"
with open(outfile, "w") as f:
    df_filtered.to_csv(f, sep=",", index=False, quoting=csv.QUOTE_ALL)
