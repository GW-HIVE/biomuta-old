import csv
import pandas as pd

def count_string_occurrences(csv_file, column_name, search_string):
    # Load the CSV file into a DataFrame
    df = pd.read_csv(csv_file)
    
    # Check if the column exists
    if column_name not in df.columns:
        print(f"Column '{column_name}' not found in the CSV file.")
        return
    
    # Count occurrences of the search string in the specified column
    count = df[column_name].str.count(search_string, na='nan').sum()

    print(f"The string '{search_string}' appears {count} times in the column '{column_name}'.")
'''
# Example usage
csv_file = '/data/shared/biomuta/generated/datasets/2024_10_22/mapping_ids/chr_pos_to_ensp_old.csv'  # Replace with your actual CSV file path
column_name = 'prot_change'  # Replace with the actual column name
search_string = 'foo'  # Replace with the string you're searching for

count_string_occurrences(csv_file, column_name, search_string)
'''


def count_unique_uniprotkb(file_path):
    """
    Counts unique values in the 'uniprotkb_canonical_ac' column up to the '-' character.

    Args:
        file_path (str): Path to the CSV file.

    Returns:
        int: Count of unique values.
    """
    unique_values = set()
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            if 'uniprotkb_canonical_ac' in row:
                value = row['uniprotkb_canonical_ac'].split('-')[0]  # Take part before "-"
                unique_values.add(value)
    return len(unique_values)

# Example usage:
file_path = "/data/shared/biomuta/downloads/biomuta.csv"
print(f"Unique uniprotkb_canonical_ac values: {count_unique_uniprotkb(file_path)}")
