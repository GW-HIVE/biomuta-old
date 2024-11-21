import pandas as pd

# Load the CSV files into DataFrames
bm = pd.read_csv('/data/shared/biomuta/generated/stats/bm_egfr_unique_combinations_sorted.csv')
cbio = pd.read_csv('/data/shared/biomuta/generated/stats/cbio_egfr_sub_reshaped.csv')
gg = pd.read_csv('/data/shared/biomuta/generated/stats/gg_egfr_extracted_columns.csv')

# Find the common rows across all three files
common_rows = pd.merge(bm, cbio, how='inner').merge(gg, how='inner')

# Find common rows between each pair of files
common_bm_cbio = pd.merge(bm, cbio, how='inner')
common_cbio_gg = pd.merge(cbio, gg, how='inner')
common_bm_gg = pd.merge(bm, gg, how='inner')

# Find the unique rows in each file
unique_bm = pd.concat([bm, common_rows, common_bm_cbio, common_bm_gg]).drop_duplicates(keep=False)
unique_cbio = pd.concat([cbio, common_rows, common_bm_cbio, common_cbio_gg]).drop_duplicates(keep=False)
unique_gg = pd.concat([gg, common_rows, common_cbio_gg, common_bm_gg]).drop_duplicates(keep=False)

# Output the results
common_rows.to_csv('/data/shared/biomuta/generated/stats/egfr_common_rows_bm_cbio_gg.csv', index=False)
unique_bm.to_csv('/data/shared/biomuta/generated/stats/unique_bm.csv', index=False)
unique_cbio.to_csv('/data/shared/biomuta/generated/stats/unique_cbio.csv', index=False)
unique_gg.to_csv('/data/shared/biomuta/generated/stats/unique_gg.csv', index=False)
common_bm_cbio.to_csv('/data/shared/biomuta/generated/stats/common_bm_cbio.csv', index=False)
common_cbio_gg.to_csv('/data/shared/biomuta/generated/stats/common_cbio_gg.csv', index=False)
common_bm_gg.to_csv('/data/shared/biomuta/generated/stats/common_bm_gg.csv', index=False)

print("Comparison complete! Check the output CSV files for results.")
