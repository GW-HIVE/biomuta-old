#!/bin/bash

# Load config.json
CONFIG_FILE="/data/shared/repos/biomuta-old/pipeline/config.json"
# Get paths from config.json
DOWNLOADS=$(jq -r '.relevant_paths.downloads' "$CONFIG_FILE")
GENERATED=$(jq -r '.relevant_paths.generated_datasets' "$CONFIG_FILE")
# Find the latest dump directory within the cbioportal directory
LATEST_DUMP=$(ls -t "$DOWNLOADS/cbioportal" | head -n1)
# Construct the output directory path using the latest dump directory name
OUTPUT_DIR="${GENERATED}/${LATEST_DUMP}/cancer_types"
# Check if the output directory already exists
if [ -d "${OUTPUT_DIR}" ]; then
    echo "Directory ${OUTPUT_DIR} already exists."
    # Ask user for confirmation to overwrite or create a new directory
    read -p "Do you want to overwrite it? (y/n): " confirm
    if [[ "$confirm" == [yY] ]]; then
        echo "Overwriting ${OUTPUT_DIR}..."
        rm -rf "${OUTPUT_DIR}"  # Delete the existing directory
        mkdir -p "${OUTPUT_DIR}"  # Recreate the directory
    else
        # Create a new directory by appending a timestamp to avoid overwriting
        TIMESTAMP=$(date +"%Y%m%d%H%M%S")
        NEW_OUTPUT_DIR="${OUTPUT_DIR}_${TIMESTAMP}"
        echo "Creating new directory: ${NEW_OUTPUT_DIR}"
        mkdir -p "${NEW_OUTPUT_DIR}"
        OUTPUT_DIR="${NEW_OUTPUT_DIR}"  # Update OUTPUT_DIR to point to the new directory
    fi
else
    # Create the output directory if it doesn't exist
    mkdir -p "${OUTPUT_DIR}"
fi
echo "Final output directory: ${OUTPUT_DIR}"