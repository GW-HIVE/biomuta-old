import pandas as pd

# Load the CSV files
cbio_table = pd.read_csv("/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/final_table.csv")
biomuta = pd.read_csv("/data/shared/biomuta/generated/stats/full_tables/biomuta_canonical.csv")

# Select the columns to compare
cbio_tuples = set(cbio_table[["uniprotkb_canonical_ac", "aa_pos", "ref_aa", "alt_aa"]].dropna().itertuples(index=False, name=None))
biomuta_tuples = set(biomuta[["uniprotkb_canonical_ac", "aa_pos", "ref_aa", "alt_aa"]].dropna().itertuples(index=False, name=None))

# Perform set operations
in_both = cbio_tuples & biomuta_tuples
only_in_cbio = cbio_tuples - biomuta_tuples
only_in_biomuta = biomuta_tuples - cbio_tuples

# Convert results to DataFrames and retrieve "source" values
in_both_df = pd.DataFrame(in_both, columns=["uniprotkb_canonical_ac", "aa_pos", "ref_aa", "alt_aa"])
in_both_df["comparison_source"] = "in_both"

only_in_cbio_df = pd.DataFrame(only_in_cbio, columns=["uniprotkb_canonical_ac", "aa_pos", "ref_aa", "alt_aa"])
only_in_cbio_df["comparison_source"] = "only_in_cbio"

only_in_biomuta_df = pd.DataFrame(only_in_biomuta, columns=["uniprotkb_canonical_ac", "aa_pos", "ref_aa", "alt_aa"])
only_in_biomuta_df["comparison_source"] = "only_in_biomuta"

# Combine results
results_df = pd.concat([in_both_df, only_in_cbio_df, only_in_biomuta_df]).sort_values(by=["Source"] + ["uniprotkb_canonical_ac", "aa_pos", "ref_aa", "alt_aa"])

# Save to a TSV file
results_df.to_csv("/data/shared/biomuta/generated/stats/full_tables/venn_4_components.tsv", sep="\t", index=False)
