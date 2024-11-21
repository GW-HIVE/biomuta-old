#!/bin/bash

# Define paths
file="/data/shared/biomuta/generated/datasets/2024_10_22/cbio_egfr_aa_changes.txt"
tmp="/data/shared/biomuta/generated/datasets/2024_10_22/cbio_egfr_tmp.txt"
subs="/data/shared/biomuta/generated/datasets/2024_10_22/cbio_egfr_sub.txt"
nonsubs="/data/shared/biomuta/generated/datasets/2024_10_22/cbio_egfr_nonsub.txt"

while IFS= read -r line
do
  # Remove both double and single quotation marks
  modified_line="${line//\"/}"
  modified_line="${modified_line//\'/}"
  
  # Write the modified line to the temporary file
  echo "$modified_line" >> "$tmp"
done < "$file"

# AA changes that are not substitutions
grep -e "fs" -e "del" -e "ins" -e "dup" -e "splice" -e "MUTATED" $tmp > $nonsubs

# AA changes that are substitutions
grep -v -e "fs" -e "del" -e "ins" -e "dup" -e "splice" -e "MUTATED" $tmp > $subs

total_n=$(wc -l < "$tmp")
nonsub_n=$(wc -l < "$nonsubs")
sub_n=$(wc -l < "$subs")
echo "Total number of lines: ${total_n}"
echo "Number of non-substitutions: ${nonsub_n}"
echo "Number of substitutions: ${sub_n}"
echo "Expected total: $(($nonsub_n + $sub_n))"

if (( total_n == nonsub_n + sub_n )); then
    echo "Count check passed."
else
    echo "Warning: Total does not match sum of substitutions and non-substitutions."
fi

rm "$tmp"