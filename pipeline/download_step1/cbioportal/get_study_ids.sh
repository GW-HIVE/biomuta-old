#!/bin/bash

# Create downloads dir: done, success
# Use API to download: done, success
# Extract whatever you need
#   name
#   studyId
#   cancerTypeId
# Testing...

# Source the utils script
THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$THIS_DIR/../../utils/utils.sh"

# Create downloads directory with today's date as name
mkdir_today

# Fetch study IDs and their respective cancer types
curl -G "https://www.cbioportal.org/api/studies" \
  -H "accept: application/json" \
  -o "${DL_TODAY}/study_ids_raw.json"

# Check if curl command was successful
if [ $? -eq 0 ]; then
    echo "Study IDs fetched successfully and saved to ${DL_TODAY}/study_ids.json"
else
    echo "Failed to fetch study IDs. Please check your internet connection or the API URL."
    exit 1
fi

# Extract fields
jq 