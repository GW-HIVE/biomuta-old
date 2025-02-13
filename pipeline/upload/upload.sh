#!/bin/bash

# Directory containing the biomuta files
input_dir="/data/shared/repos/biomuta-old/generated_datasets/compiled"

# Base filename pattern
base_filename="biomuta_v"

# Upload directory
output_base_dir="/data/shared/repos/biomuta-old/nginx-file-server/data"

# Find the most recent file
latest_file=$(ls -1v "$input_dir"/"$base_filename"*.csv | tail -n 1)

# Check if a file was found
if [[ -z "$latest_file" ]]; then
    echo "No files matching the pattern '$base_filename*.csv' found in $input_dir."
    exit 1
fi

# Extract the version number from the filename
version=$(basename "$latest_file" | grep -oP '(?<=_v)\d+(?=\.csv)')

# Create the new directory for the version
new_dir="$output_base_dir/$version.0"
mkdir -p "$new_dir"

# Copy the latest file to the new directory
cp "$latest_file" "$new_dir/biomuta.csv"

# Update the "current" symlink to point to the new directory
ln -sfn "$new_dir" "$output_base_dir/current"

echo "Latest file: $latest_file"
echo "Copied to: $new_dir"
echo "Symlink 'current' updated to point to: $new_dir"
