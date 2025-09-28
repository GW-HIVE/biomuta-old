#!/bin/bash

# Create the truncated directory if it doesn't exist
truncated_dir="/data/shared/biomuta/downloads/cbioportal/current/mutations/truncated/"
mkdir -p "$truncated_dir"

for file in /data/shared/biomuta/downloads/cbioportal/current/mutations/*.json; do
    if ! tail -c 10 "$file" | grep -q '}'; then
        echo "Potentially truncated: $file"
        mv "$file" "$truncated_dir"
    fi
done