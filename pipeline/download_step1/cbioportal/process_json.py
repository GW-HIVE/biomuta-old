import json
import os
import pandas as pd

def get_headers_from_json(file_path):
    """
    Extract the column names from a JSON file.
    """
    try:
        with open(file_path, 'r') as file:
            data = json.load(file)
            headers = set(data[0].keys()) if data else set()
        return data, headers
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON from file {file_path}: {e}")
        return [], set()  # Return an empty list and set if there's an error
    

def find_common_and_unique_headers(json_files):
    all_data = []
    file_headers = {}

    # Collect data and headers from each file
    for file in json_files:
        data, headers = get_headers_from_json(file)
        if data:  # Proceed only if data was loaded successfully
            all_data.append(data)
            file_headers[file] = headers

    # Find common headers across all files
    common_headers = set.intersection(*[h for h in file_headers.values()])

    # Find unique headers for each file
    unique_headers = {file: headers - common_headers for file, headers in file_headers.items()}

    return common_headers, unique_headers, all_data

def create_excel_files(common_headers, unique_headers, all_data, output_dir="."):
    # Create a DataFrame for common headers
    common_data_frames = []
    
    for data in all_data:
        if data:
            df = pd.DataFrame(data)
            common_data_frames.append(df[common_headers])

    # Concatenate all common data frames and save to Excel
    common_df = pd.concat(common_data_frames, ignore_index=True)
    common_df.to_excel(os.path.join(output_dir, "common_headers_data.xlsx"), index=False)

    print(f"Excel file with common headers created: {output_dir}/common_headers_data.xlsx")

    # Create Excel files for each JSON file with unique headers
    for file, headers in unique_headers.items():
        if headers:
            # Create a DataFrame for the unique headers
            df = pd.DataFrame(get_headers_from_json(file)[0])  # Load data
            unique_df = df[list(headers)]
            unique_file_name = os.path.basename(file).replace('.json', '_unique_headers.xlsx')
            unique_df.to_excel(os.path.join(output_dir, unique_file_name), index=False)
            print(f"Excel file created for unique headers: {output_dir}/{unique_file_name}")

# Example usage
if __name__ == "__main__":
    # Folder containing your JSON files
    folder_path = "/data/shared/pipelines/cbioportal/mutations"
    output_dir = "/data/shared/pipelines/cbioportal/excel"  # Specify your output directory here
    
    # List all json files in the folder
    json_files = [os.path.join(folder_path, file) for file in os.listdir(folder_path) if file.endswith('.json')]

    # Find common and unique headers
    common_headers, unique_headers, all_data = find_common_and_unique_headers(json_files)

    # Create Excel files
    create_excel_files(common_headers, unique_headers, all_data, output_dir)
