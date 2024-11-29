import os
import ijson

def find_incomplete_json_files(directory, output_file="incomplete_files.txt"):
    # List to store the names of incomplete files
    incomplete_files = []

    # Iterate over all files in the directory
    for filename in os.listdir(directory):
        if filename.endswith(".json"):  # Process only JSON files
            filepath = os.path.join(directory, filename)
            try:
                # Try parsing the JSON file incrementally
                with open(filepath, "r") as f:
                    for _ in ijson.items(f, ""):
                        pass  # If it parses successfully, do nothing
            except (ijson.JSONError, ValueError):
                # If parsing fails, the file is incomplete or corrupted
                incomplete_files.append(filename)

    # Write the list of incomplete files to the output file
    with open(output_file, "w") as f:
        for filename in incomplete_files:
            f.write(filename + "\n")

    print(f"Found {len(incomplete_files)} incomplete files. Results saved to '{output_file}'.")

# Specify the directory containing your JSON files
directory_path = "/data/shared/biomuta/downloads/cbioportal/2024_10_21/mutationss"

# Run the function
find_incomplete_json_files(directory_path)
