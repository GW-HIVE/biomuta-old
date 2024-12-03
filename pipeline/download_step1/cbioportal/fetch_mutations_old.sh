#!/bin/bash

# Load paths from config.json
CONFIG_FILE="/path/to/config.json"
DOWNLOADS_DIR=$(jq -r '.relevant_paths.generated_datasets' "$CONFIG_FILE")

# Today's date
TODAY=$(date +"%Y_%m_%d")

# Output directory
OUTPUT_DIR="${DOWNLOADS_DIR}/${TODAY}"
mkdir -p "${OUTPUT_DIR}/mutations"

# Base URLs for the API
STUDY_URL="https://www.cbioportal.org/api/studies"
MUTATIONS_URL="https://www.cbioportal.org/api/molecular-profiles"

# Get study IDs
curl -G "${STUDY_URL}" \
     -H "accept: application/json" \
     -o "${OUTPUT_DIR}/all_studies.json"

#original command curl -G "https://www.cbioportal.org/api/studies" -H "accept: application/json" -o all_studies.json

#extract study IDs
jq -r '.[].studyId' "${OUTPUT_DIR}/all_studies.json" > "${OUTPUT_DIR}/study_ids.txt"

# File containing study IDs, one per line
STUDY_IDS_FILE="${OUTPUT_DIR}/study_ids.txt"

# Mutations endpoint requires Molecular Profile IDs and Sample List IDs

# Iterate over each study ID in the file
while IFS= read -r study_id; do
  echo "Processing study ID: $study_id"

  # Fetch molecular profile IDs for the current study
  curl -G "${STUDY_URL}/${study_id}/molecular-profiles" \
       -H "accept: application/json" \
       -o "${OUTPUT_DIR}/${study_id}_molecular_profiles.json"

  if [ $? -eq 0 ]; then
    echo "  Molecular profiles fetched successfully for study ID: $study_id"
  else
    echo "  Failed to fetch molecular profiles for study ID: $study_id"
    continue
  fi

  # Fetch sample list IDs for the current study
  curl -G "${STUDY_URL}/${study_id}/sample-lists" \
       -H "accept: application/json" \
       -o "${OUTPUT_DIR}/${study_id}_sample_lists.json"

  if [ $? -eq 0 ]; then
    echo "  Sample lists fetched successfully for study ID: $study_id"
  else
    echo "  Failed to fetch sample lists for study ID: $study_id"
    continue
  fi

  # Extract molecularProfileId from the JSON response
  molecular_profile_ids=$(jq -r '.[].molecularProfileId' "${OUTPUT_DIR}/${study_id}_molecular_profiles.json")
  # Extract sampleListId from the JSON response
  sample_list_ids=$(jq -r '.[].sampleListId' "${OUTPUT_DIR}/${study_id}_sample_lists.json")

  # Iterate over each molecular profile ID and sample list ID
  for molecular_profile_id in $molecular_profile_ids; do
    for sample_list_id in $sample_list_ids; do
      echo "  Fetching mutations for molecular profile ID: $molecular_profile_id and sample list ID: $sample_list_id"

      # Fetch mutations for the molecular profile ID and sample list ID
      curl -G "${MUTATIONS_URL}/${molecular_profile_id}/mutations" \
           -H "accept: application/json" \
           --data-urlencode "sampleListId=${sample_list_id}" \
           -o "${OUTPUT_DIR}/mutations/${molecular_profile_id}_${sample_list_id}.json"

      # Check if mutation data fetch was successful
      if [ $? -eq 0 ]; then
        echo "    Mutations fetched successfully for $molecular_profile_id with sample list $sample_list_id"
      else
        echo "    Failed to fetch mutations for $molecular_profile_id with sample list $sample_list_id"
      fi

      # Add a small delay to avoid hitting rate limits
      sleep 5
    done
  done

done < "$STUDY_IDS_FILE"
