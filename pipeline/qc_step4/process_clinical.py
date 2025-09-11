import csv
import pandas as pd

def process_csv_files(clinical_file, biomuta_file):
    """
    Process the CSV files according to the requirements:
    1. Match rows using matching headers, removing unmatched rows from clinical_file
    2. Keep only 3 columns and rename as specified
    3. Split into 2 files
    """
    print("Loading CSV files...")
    clinical_df = pd.read_csv(clinical_file)
    biomuta_df = pd.read_csv(biomuta_file)

    print(f"Original clinical data shape: {clinical_df.shape}")
    print(f"Original biomuta data shape: {biomuta_df.shape}")

    # Step 1: match rows using common headers (chr_id, ref_nt, alt_nt, aa_pos)
    common_columns = ['chr_id', 'ref_nt', 'alt_nt', 'aa_pos']
    biomuta_df["_match_id"] = biomuta_df.groupby(common_columns).cumcount() # Assign a unique index to each duplicate row in biomuta_df within each duplicate group
    clinical_df["_match_id"] = clinical_df.groupby(common_columns).cumcount() # Assign a correspoding unique index to each duplicate row in clinical_df within each duplicate group
    print("Matching rows based on common columns with one-to-one matching...")
    merged_df = pd.merge( # Create a merged dataframe by matching rows
        clinical_df,
        biomuta_df,
        on=common_columns + ["_match_id"], # One-to-one matching
        how='inner'
    )
    print(f"After merging, clinical data shape: {merged_df.shape}")

    # Step 2: rename columns
    ## Create temporary DataFrame with correct column names for splitting
    temp_sample_df = merged_df[['sample_id', 'proj__project_id']].copy()
    temp_sample_df.rename(columns={'sample_id': 'submitter_id'}, inplace=True)
    print(f"sample_df has {len(temp_sample_df)} rows")
    print("sample_df first 10 rows:")
    print(temp_sample_df.head(10))

    temp_patient_df = merged_df[['patient_id', 'proj__project_id']].copy()
    temp_patient_df.rename(columns={'patient_id': 'submitter_id'}, inplace=True)
    print(f"patient_df has {len(temp_sample_df)} rows")
    print("patient_df first 10 rows:")
    print(temp_patient_df.head(10))

    # Step 3: save the two files
    temp_sample_df.to_csv('/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/clinical-information-sample.csv', index=False, quoting=csv.QUOTE_ALL)
    temp_patient_df.to_csv('/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/clinical-information-patient.csv', index=False, quoting=csv.QUOTE_ALL)

# Usage
process_csv_files("/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/clinical-information.csv", "/data/shared/repos/biomuta-old/generated_datasets/compiled/biomuta_v6.1.csv")