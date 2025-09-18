# 2_parse_gff.py Documentation (auto-generated, unreviewed)

## Overview

This script maps chromosomal positions from BED files to Ensembl protein IDs (ENSP) using GFF3 annotation data. It serves as a key component in a bioinformatics pipeline for cancer variant mapping, converting genomic coordinates to protein identifiers for downstream analysis.

## Purpose

The script processes genomic variant data by:
- Loading pre-built GFF3 database containing Ensembl annotations
- Mapping chromosomal positions to protein-coding sequences (CDS)
- Extracting ENSP (Ensembl Protein) IDs for variant positions
- Converting BED format input to TSV output with protein mappings

## Dependencies

```python
import csv
import gffutils
import logging
import sys
from pathlib import Path
```

### External Dependencies
- **gffutils**: Python package for working with GFF/GTF files
- **Custom utilities**: Local utils module for configuration and paths

## Configuration

### Logging Setup
- **Log file**: `cancer_mapping.log`
- **Format**: Timestamp, log level, and message
- **Level**: INFO and above
- **Mode**: Append mode for persistent logging

### File Paths
- **GFF database**: `/data/shared/repos/biomuta-old/downloads/ensembl/Homo_sapiens.GRCh38.113.db`
- **Input BED file**: Retrieved from configuration via `hg38_combined.bed`
- **Output file**: `chr_pos_to_ensp.tsv` in mapping_ids directory

## Core Functions

### `get_ensp_for_position(chrom, start, end)`

Maps chromosomal coordinates to Ensembl protein IDs.

**Parameters:**
- `chrom` (str): Chromosome identifier (e.g., 'chr1', 'chrX')
- `start` (int): Start position (0-based)
- `end` (int): End position (0-based)

**Returns:**
- `list`: ENSP IDs found at the position, or ['N/A'] if none found

**Process:**
1. Normalizes chromosome format (removes 'chr' prefix)
2. Converts to 1-based coordinates for GFF compatibility
3. Queries database for CDS features in the specified region
4. Extracts protein_id attributes from matching features

### `clean_chr_id(chr_id)`

Standardizes chromosome identifiers for consistent output.

**Parameters:**
- `chr_id` (str): Raw chromosome identifier

**Returns:**
- `str`: Cleaned chromosome ID

**Mapping Rules:**
- Removes 'chr' prefix
- Maps 'X' → '23'
- Maps 'Y' → '24'
- Preserves numeric chromosomes as-is

### `process_bed_file(input_file, output_file)`

Main processing function that converts BED file to ENSP-mapped TSV.

**Parameters:**
- `input_file` (str): Path to input BED file
- `output_file` (str): Path for output TSV file

**Process:**
1. Reads tab-separated BED file
2. Skips headers and empty rows
3. For each variant position:
   - Extracts genomic coordinates
   - Queries for ENSP IDs
   - Writes results with one row per ENSP ID

## Input Format (BED File)

Expected tab-separated columns:
```
chr_id    start_pos    end_pos    entrez_gene_id    prot_change
```

**Example:**
```
chr10     43163970     43163971   12345            p.Val123Met
```

## Output Format (TSV File)

Generated columns:
```
chr_id    start_pos    end_pos    entrez_gene_id    prot_change    ensp
```

**Example:**
```
10        43163970     43163971   12345            p.Val123Met    ENSP00000123456
```

## Database Structure

### GFF3 Database Creation (Commented Code)
The script includes commented code for initial database creation:

```python
# One-time database creation from GFF3 file
db = gffutils.create_db(
    gff_file,
    database,
    id_spec=['ID', 'Name', 'Parent'],
    merge_strategy='create_unique',
    force=True,
    keep_order=True
)
```

**Parameters Explained:**
- `id_spec`: Hierarchy for feature ID assignment
- `merge_strategy`: Handles duplicate IDs by creating unique suffixes
- `force`: Overwrites existing database
- `keep_order`: Maintains original GFF3 feature order

### Feature Types
The database contains various genomic features, with the script specifically querying 'CDS' (Coding Sequence) features to find protein-coding regions.

## Error Handling

- **Missing Features**: Returns ['N/A'] when no ENSP IDs found
- **Empty Rows**: Skips empty or header rows in input
- **Data Validation**: Converts positions to integers for proper querying

## Performance Considerations

- **Database Loading**: Uses pre-built SQLite database for fast queries
- **Memory Efficiency**: Processes files line-by-line rather than loading entirely
- **Query Optimization**: Targets specific feature types (CDS) for relevant results

## Usage Example

The script is configured to run with specific input/output paths:

```python
config_obj = get_config()
wd = Path(config_obj["relevant_paths"]["generated_datasets"])
input_file = wd / '2024_10_22' / 'liftover' / 'hg38_combined.bed'
output_file = wd / '2024_10_22' / 'mapping_ids' / 'chr_pos_to_ensp.tsv'

process_bed_file(input_file, output_file)
```

## Integration Notes

This script appears to be part of a larger bioinformatics pipeline:
- **Step 1**: Likely genome coordinate liftover (based on input path)
- **Step 2**: This script - genomic position to protein mapping
- **Step 3+**: Downstream analysis using ENSP identifiers

## Genome Build

- **Reference**: Human genome GRCh38 (Ensembl release 113)
- **Coordinate System**: 0-based input (BED format), converted to 1-based for GFF queries

## Logging Output

The script logs important events and results:
- Database loading status
- Query results for debugging
- Processing progress information

## Maintenance Notes

- **Database Updates**: Requires periodic updates with new Ensembl releases
- **Path Configuration**: Uses configuration system for flexible file paths
- **Feature Types**: Currently hardcoded to query 'CDS' features only