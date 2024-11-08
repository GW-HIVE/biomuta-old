#!/bin/bash

# This script extracts study IDs and the corresponding cancer names. It takes as input the output of cancer_types.sh

# Define the directory where your JSON files are located
input_dir="/data/shared/biomuta/downloads/cbioportal/2024_10_21/cancer_types"

# Define the output TSV file
output_file="/data/shared/biomuta/generated/datasets/2024_10_22/cancer_type_per_study.json"

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

# Make a list of unique cancer names in json format
jq -r '.[].cancerType' $output_file | sort | uniq | jq -R . | jq -s . > unique_cancer_names.json