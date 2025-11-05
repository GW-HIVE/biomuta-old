#!/bin/bash

# Get the directory of this script
UTILS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="$UTILS_DIR/../pipeline/config.json"

# Export paths from the config file
export DOWNLOADS=$(jq -r '.relevant_paths.downloads' "$CONFIG_FILE")
export GENERATED=$(jq -r '.relevant_paths.generated_datasets' "$CONFIG_FILE")

# Function to create a directory with today's date as name
mkdir_today() {
    TODAY=$(date +"%Y_%m_%d") # Get today's date in yyyy_mm_dd format
    DL_TODAY="${DOWNLOADS}/cbioportal/${TODAY}"
    MUTATIONS_DIR="${DL_TODAY}/mutations"

    # Check if today's download directory already exists
    if [ -d "${DL_TODAY}" ]; then
        echo "Directory ${DL_TODAY} already exists."
        # Ask user for confirmation to overwrite
        read -p "Do you want to overwrite it? (y/n): " confirm
        if [[ "$confirm" == [yY] ]]; then
            echo "Overwriting ${DL_TODAY}..."
            rm -rf "${DL_TODAY}"  # Delete the existing directory
            mkdir -p "${DL_TODAY}" "${MUTATIONS_DIR}"  # Recreate the directory
            echo "Directories created: ${DL_TODAY}, ${MUTATIONS_DIR}"
        elif [[ "$confirm" == [nN] ]]; then
            echo "Exiting without making changes."
            exit 0
        else
            echo "Invalid choice. Exiting."
            exit 0
        fi
    else
        mkdir -p "${DL_TODAY}" "${MUTATIONS_DIR}"
        echo "Directories created: ${DL_TODAY}, ${MUTATIONS_DIR}"
    fi
}

# Function to get the latest directory in a path
get_latest_dir() {
    local base_path=$1
    ls -t "$base_path" | tail -n1
}

# Find the latest dump cbioportal directory
LATEST_DUMP=$(get_latest_dir "$DOWNLOADS/cbioportal")
# Latest generated data directory
LATEST_GEN=$(get_latest_dir "$GENERATED")