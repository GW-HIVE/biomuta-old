#!/bin/bash

bm='/data/shared/biomuta/generated/stats/bm_egfr_unique_combinations_sorted.csv'
cbio='/data/shared/biomuta/generated/stats/cbio_egfr_sub_reshaped.csv'
out='/data/shared/biomuta/generated/stats/bm_cbio_egfr_diff_full.txt'
unique_to_bm='/data/shared/biomuta/generated/stats/unique_to_bm.csv'

sorted_bm=$(sort <(sed 's/^[[:space:]]*//;s/[[:space:]]*$//' $bm))
sorted_cbio=$(sort <(sed 's/^[[:space:]]*//;s/[[:space:]]*$//' $cbio))

# Create temporary files
tmp_bm=$(mktemp)
tmp_cbio=$(mktemp)

# Write sorted contents to the temporary files
echo "$sorted_bm" > "$tmp_bm"
echo "$sorted_cbio" > "$tmp_cbio"

# Run comm command and redirect output to $out
comm "$tmp_bm" "$tmp_cbio" > "$out"

# Clean up temporary files
rm "$tmp_bm" "$tmp_cbio"

# Count lines in the first column (unique to bm)
awk -F'\t' '{if ($1 != "") print $0}' $out > $unique_to_bm
first_column_count=$(wc -l < "$unique_to_bm")

# Count lines in the second column (unique to cbio)
second_column_count=$(awk -F'\t' '{if ($2 != "") print $0}' $out | wc -l)

# Count lines in the third column (common lines)
third_column_count=$(awk -F'\t' '{if ($3 != "") print $0}' $out | wc -l)

# Output the results
echo "Lines in the first column (unique to bm): $first_column_count"
echo "Lines in the second column (unique to cbio): $second_column_count"
echo "Lines in the third column (common lines): $third_column_count"
