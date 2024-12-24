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
    mkdir -p "${DL_TODAY}" "${MUTATIONS_DIR}"
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