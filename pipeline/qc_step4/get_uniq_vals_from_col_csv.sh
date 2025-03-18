#!/bin/bash

# Require to specify filename and column index
if [ $# -lt 2 ]; then
        echo ''
        echo "Error: Please specify the filename and column index."
        echo "Usage: $0 <filename> <column_index>"
        echo ''
        exit 1
fi

# Initialize filename and column index
file=$1
col=$2

# Define paths
OUTFILE="/data/shared/biomuta/generated/qc/uniq_vals_col_$col.csv"

awk -F',' -v col="$col" 'FNR==NR{a[$col]++;next} a[$col]==1' "$file" "$file" > "$OUTFILE"

echo "Unique values written to: $OUTFILE"