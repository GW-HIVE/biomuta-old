import pandas as pd

def open_csv_as_df(file_path):
    df = pd.read_csv(file_path)
    return df

def increment_column_values(df, column_name):
    df[column_name] = df[column_name] + 1
    return df

def write_df_to_csv(df, file_path, index=False):
    df.to_csv(file_path, index=index)

df = open_csv_as_df('/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/mapping_ids/filtered_chr_pos_with_uniprot.csv')
incremented_df = increment_column_values(df, 'start_pos')
write_df_to_csv(incremented_df, '/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/mapping_ids/1-based_chr_pos.csv')