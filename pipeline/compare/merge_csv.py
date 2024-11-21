import pandas as pd

# Load both CSV files
file = '/data/shared/repos/biomuta-old/downloads/glygen/snv_disease_associated_mutations_P00533-1.csv'
df = pd.read_csv(file)

# Specify the column and value to filter out
target_column = 'Source Xref DB'
specific_value = 'BioMuta'

# Remove rows where the target column has the specific value
result_df = df[df[target_column] != specific_value]

# Specify the columns to consider for identifying duplicates
duplicate_columns = ['Start Position', 'End Pos', 'Sequence']

# Remove duplicates based on the specified columns
result_df = result_df.drop_duplicates(subset=duplicate_columns)

def extract_columns(arg_df, output_file):
    # Select only the 2nd and 3rd columns
    df_two_three = arg_df.iloc[:, 1:3]
    # Split the 3rd column and create two new columns
    df_two_three[['ref_aa', 'alt_aa']] = df_two_three['Sequence'].str.split(' -> ', expand=True)
    # Drop the original column if no longer needed
    df_two_three = df_two_three.drop(columns=['Sequence'])
    # Save to a new CSV file
    df_two_three.to_csv(output_file, index=False)

output_file = '/data/shared/biomuta/generated/stats/gg_egfr_extracted_columns.csv'
extract_columns(result_df, output_file)

print("Output saved.")
