#!/bin/bash


# Get this script's directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Load config.json
CONFIG_FILE="$SCRIPT_DIR/../../config.json"


# Get paths from config
input_dir=$(jq -r '.relevant_paths.downloads + "/cbioportal/2024_10_21/cancer_types"' "$CONFIG_FILE")
output_dir=$(jq -r '.relevant_paths.generated_datasets + "/2024_10_22"' "$CONFIG_FILE")

# Define the output files
output_file="$output_dir/cancer_type_per_study.json"
unique_cancer_names_file="$output_dir/unique_cancer_names.json"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Initialize the JSON array
echo "[" > "$output_file"

# Track if this is the first record (for proper JSON formatting)
first_record=true

# Iterate over each JSON file in the directory
for file in "$input_dir"/*.json; do
  # Extract studyId and cancerType.name using jq
  study_id=$(jq -r '.studyId' "$file")
  cancer_type=$(jq -r '.cancerType.name' "$file")

  # Format the data as a JSON object
  json_object=$(jq -n \
    --arg sid "$study_id" \
    --arg ctype "$cancer_type" \
    '{studyId: $sid, cancerType: $ctype}')

  # Append a comma if it's not the first record
  if [ "$first_record" = true ]; then
    first_record=false
  else
    echo "," >> "$output_file"
  fi

  # Append the JSON object to the file
  echo "$json_object" >> "$output_file"
done

# Close the JSON array
echo "]" >> "$output_file"

echo "Data successfully written to $output_file"

# Make a list of unique cancer names in JSON format
jq -r '.[].cancerType' "$output_file" | sort | uniq | jq -R . | jq -s . > "$unique_cancer_names_file"

echo "Unique cancer names written to $unique_cancer_names_file"
