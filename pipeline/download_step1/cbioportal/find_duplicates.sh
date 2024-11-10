#!/bin/bash
awk '{print $2}' duplicate_files.txt | grep -v '^$' | while read -r file; do
    # Ensure the file exists before calculating its checksum
    if [[ -f "$file" ]]; then
        checksum=$(md5sum "$file" | cut -d' ' -f1)
        
        # Find all files with the same checksum
        duplicate_files=($(grep "^$checksum" duplicate_files.txt | awk '{print $2}'))
        
        # Check if we have found duplicates
        if [[ ${#duplicate_files[@]} -gt 0 ]]; then
            first_file=${duplicate_files[0]}
            
            # Delete the file if it's not the first occurrence
            if [[ "$file" != "$first_file" ]]; then
                echo "Would delete: $file"
                # Uncomment the next line to actually delete the file
                rm "$file"
            fi
        else
            echo "No duplicates found for: $file"
        fi
    else
        echo "File does not exist: $file"
    fi
done