# 1_generate_cancer_do_json.py (auto-generated, unreviewed)
Cancer Disease Ontology (DO) Mapping Tool

A Python script for mapping cancer types to Disease Ontology terms using hierarchical mapping sources.

## Overview

This tool processes cancer type data and maps it to standardized Disease Ontology (DO) terms. It uses a two-tier mapping approach with primary sources (CIVIC and COSMIC) and fallback mappings to ensure comprehensive coverage.

## Features

- **Hierarchical Mapping**: Primary mapping from CIVIC/COSMIC sources with fallback support
- **Interactive Confirmation**: User-prompted directory selection for safety
- **Comprehensive Logging**: Detailed logging with configurable levels
- **JSON Processing**: Handles multiple JSON input/output formats
- **Study-Level Mapping**: Maps individual studies to DO terms via cancer types

## Prerequisites

- Python 3.7+
- Required files:
  - `config.json` (configuration file with paths)
  - `combined_do_mapping.json` (primary CIVIC/COSMIC mappings)
  - `fallback_cbio_doid_mapping.json` (fallback mappings)
  - `unique_cancer_names.json` (input cancer types)
  - `cancer_type_per_study.json` (study-cancer type associations)

## Configuration

The script expects a `config.json` file with the following structure:

```json
{
  "relevant_paths": {
    "mapping": "/path/to/mapping/files",
    "generated_datasets": "/path/to/datasets"
  }
}
```

## Input Files

### `unique_cancer_names.json`
```json
[
  "breast cancer",
  "lung adenocarcinoma",
  "glioblastoma"
]
```

### `cancer_type_per_study.json`
```json
[
  {
    "studyId": "study001",
    "cancerType": "breast cancer"
  }
]
```

### `combined_do_mapping.json`
```json
{
  "civic": [
    {
      "cancer_name": "breast cancer",
      "do_term": "DOID:1612"
    }
  ],
  "cosmic": [
    {
      "cancer_name": "lung adenocarcinoma",
      "do_term": "DOID:3910"
    }
  ]
}
```

## Output Files

### `cancer_types_with_do.json`
Maps each cancer type to its corresponding DO term:
```json
[
  {
    "cancerType": "breast cancer",
    "do_name": "DOID:1612"
  }
]
```

### `study_ids_with_do.json`
Maps each study ID to its corresponding DO term:
```json
[
  {
    "studyId": "study001",
    "do_name": "DOID:1612"
  }
]
```

## Usage

1. **Prepare configuration**: Ensure `config.json` is properly configured
2. **Place input files**: Put required JSON files in the appropriate directories
3. **Run the script**:
   ```bash
   python 1_generate_cancer_do_json.py
   ```
4. **Confirm directory**: The script will prompt for confirmation of the latest dataset directory
5. **Review outputs**: Check generated JSON files and log file for results

## Mapping Logic

The script uses a hierarchical approach:

1. **Primary Mapping**: Searches CIVIC and COSMIC mappings for exact matches (case-insensitive)
2. **Fallback Mapping**: Uses keyword-based matching if no primary match found
   - Exact keyword match takes precedence
   - Word-level keyword matching as secondary option
   - Handles multiple matches with warnings

## Logging

The script generates detailed logs in `cancer_mapping.log`:

- **INFO**: Successful mappings and file operations
- **WARNING**: Multiple matches or unmapped cancer types
- **ERROR**: File access or processing errors

### Log Levels
- `ERROR`: Critical failures
- `WARNING`: Potential issues requiring attention
- `INFO`: General operations and successful mappings
- `DEBUG`: Detailed debugging information

## Error Handling

- **Multiple Matches**: Returns "Too many cancer type matches." with warning
- **No Matches**: Returns "NA" with warning log
- **Missing Files**: Script will fail with file not found errors
- **User Abort**: Clean exit when user declines directory confirmation

## File Structure

```
project/
├── config.json
├── 1_generate_cancer_do_json.py
├── cancer_mapping.log
├── mapping/
│   ├── combined_do_mapping.json
│   └── fallback_cbio_doid_mapping.json
└── generated_datasets/
    └── [timestamp_directory]/
        ├── unique_cancer_names.json
        ├── cancer_type_per_study.json
        ├── cancer_types_with_do.json
        └── study_ids_with_do.json
```

## Best Practices

- **Version Control**: Add `*.log` to `.gitignore`
- **Backup**: Keep backups of mapping files before updates
- **Validation**: Review log files after each run
- **Testing**: Test with small datasets before full processing

## Troubleshooting

### Common Issues

1. **FileNotFoundError**: Check config.json paths and file existence
2. **Permission Errors**: Ensure write permissions for output directory
3. **JSON Decode Errors**: Validate JSON file formatting
4. **Empty Results**: Check input file format and mapping file completeness

### Debug Steps

1. Check log file for detailed error messages
2. Verify all input files exist and are properly formatted
3. Confirm config.json paths are correct
4. Test with a small subset of data first

## Contributing

When contributing to this tool:

1. Maintain logging consistency
2. Add appropriate error handling
3. Update documentation for new features
4. Test with various input formats
5. Follow existing code style

## Notes

- The script automatically selects the most recently created dataset directory
- Case-insensitive matching is used for cancer type comparison
- Multiple keyword matches in fallback mapping trigger warnings
- User confirmation is required before processing begins

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