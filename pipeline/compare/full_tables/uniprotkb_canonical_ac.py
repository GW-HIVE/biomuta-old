import pandas as pd

# Load the CSV files
cbio_table = pd.read_csv("/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/final_table.csv")
biomuta = pd.read_csv("/data/shared/biomuta/downloads/biomuta_with_masterlist_check_and_isoforms.csv")

# Filter biomuta for "Found in masterlist"
biomuta_canonical = biomuta[biomuta["isoform_status"] == "Found in masterlist"]

# Extract relevant columns and remove duplicates
cbio_set = set(cbio_table["uniprotkb_canonical_ac"].dropna())
biomuta_set = set(biomuta_canonical["uniprotkb_canonical_ac"].dropna())

in_both = cbio_set & biomuta_set
only_in_cbio = cbio_set - biomuta_set
only_in_biomuta = biomuta_set - cbio_set

# Create DataFrames for each category
in_both_df = pd.DataFrame({"uniprotkb_canonical_ac": list(in_both), "Source": "in_both"})
only_in_cbio_df = pd.DataFrame({"uniprotkb_canonical_ac": list(only_in_cbio), "Source": "only_in_cbio"})
only_in_biomuta_df = pd.DataFrame({"uniprotkb_canonical_ac": list(only_in_biomuta), "Source": "only_in_biomuta"})

# Combine results
results_df = pd.concat([in_both_df, only_in_cbio_df, only_in_biomuta_df]).sort_values(by=["Source", "uniprotkb_canonical_ac"])

# Save to a TSV file
results_df.to_csv("/data/shared/biomuta/generated/stats/full_tables/uniprotkb_canonical_ac_comparison.tsv", sep="\t", index=False)
