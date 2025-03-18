import pandas as pd

def process_csv_files(clinical_file, biomuta_file):
    """
    Process the CSV files according to the requirements:
    1. Match rows using matching headers, removing unmatched rows from clinical_file
    2. Split into 2 files, keeping only 3 columns in clinical_file
    3. Rename columns as specified
    """
    print("Loading CSV files...")
    clinical_df = pd.read_csv(clinical_file)
    biomuta_df = pd.read_csv(biomuta_file)

    print(f"Original clinical data shape: {clinical_df.shape}")
    print(f"Original biomuta data shape: {biomuta_df.shape}")

    # Step 1: match rows using common headers (chr_id, ref_nt, alt_nt, aa_pos)
    common_columns = ['chr_id', 'ref_nt', 'alt_nt', 'aa_pos']
    print("Matching rows based on common columns...")
    merged_df = pd.merge( # Create a merged dataframe by matching rows
        clinical_df,
        biomuta_df[common_columns].drop_duplicates(),
        on=common_columns,
        how='inner'
    )
    print(f"After merging, clinical data shape: {merged_df.shape}")

    # Step 2: copy and only keep the 3 required columns
    final_df = merged_df[['sample_id', 'patient_id', 'proj__project_id']].copy()

    # Step 3: rename columns
    ## Create temporary DataFrame with correct column names for splitting
    temp_sample_df = merged_df[['sample_id', 'proj__project_id']].copy()
    temp_sample_df.rename(columns={'sample_id': 'submitter_id'}, inplace=True)
    temp_patient_df = merged_df[['patient_id', 'proj__project_id']].copy()
    temp_patient_df.rename(columns={'patient_id': 'submitter_id'}, inplace=True)
    ## Save the two files
    temp_sample_df.to_csv('/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/clinical-information-sample.csv', index=False)
    temp_patient_df.to_csv('/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/clinical-information-patient.csv', index=False)