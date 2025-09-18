# 1_chr_pos_to_bed.py (auto-generated, only inputs and outputs, not markdown)
Inputs

Input directory: A folder containing JSON files with mutation data from cBioPortal

Location: downloads/cbioportal/current/mutations/
Files: Multiple .json files containing genomic variant records
Each JSON file contains an array of mutation records with fields like:

ncbiBuild (genome build: GRCh37, GRCh38, etc.)
variantType (SNP, DEL, INS, etc.)
chr (chromosome)
startPosition and endPosition
entrezGeneId
proteinChange





Outputs

Output file: A single BED format file

Location: generated_datasets/current/liftover/hg19entrez_build_protChange.bed
Format: Tab-separated values with 6 columns:

Chromosome (with 'chr' prefix, X/Y converted from 23/24)
Start position (0-based, converted from 1-based)
End position


Entrez Gene ID


NCBI Build version
Protein change description





Data Filtering
The script applies several filters to select only:

SNP variants (single nucleotide polymorphisms)
Records with valid genome builds (excludes 'NA')
Non-mitochondrial chromosomes (excludes 'MT')
Non-splice site mutations
Unique records (duplicates removed)

The script processes the JSON files in batches and outputs progress information, reporting the total number of files processed and any errors encountered.

# 2_liftover.sh (auto-generated, unreviewed)

## Overview
This bash script performs genomic coordinate liftover from GRCh37 (hg19) to GRCh38 coordinate systems for biomutation data. It processes BED files containing genomic positions and protein changes, using multiple liftover strategies to maximize successful coordinate conversions.

## Purpose
The script converts genomic coordinates from the older GRCh37/hg19 reference genome to the newer GRCh38 reference genome, which is essential for maintaining consistency with current genomic databases and analysis pipelines.

## Input Requirements
- **Primary input file**: `hg19entrez_build_protChange.bed` 
  - Expected format: chromosome, start position, end position, entrez gene ID, protein change, genome build
  - Contains mixed genome builds (some already in GRCh38, others in GRCh37/hg19)
- **Chain files** (for coordinate conversion):
  - `ucscHg19ToHg38.over.chain` - UCSC liftover chain file
  - `ensembl_GRCh37_to_GRCh38.chain` - Ensembl liftover chain file
- **External tool**: `liftOver` executable must be present in the current directory

## Workflow Steps

### 1. Data Separation
```bash
# Separate GRCh38 records (already in target format)
awk '$5 == "GRCh38"' input.bed → cbio_hg38.bed

# Separate GRCh37/hg19 records (need conversion)  
awk '$5 != "GRCh38"' input.bed → hg19entrez_protChange.bed
```

### 2. Two-Stage Liftover Process
The script employs a cascading liftover strategy to maximize successful conversions:

**Stage 1: UCSC Chain File**
- Converts GRCh37 coordinates using UCSC's chain file
- Successfully mapped coordinates → `ucsc_hg38entrez_protChange.bed`
- Failed mappings → `ucsc_unmapped_entrez_protChange.bed`

**Stage 2: Ensembl Chain File**
- Processes remaining unmapped coordinates using Ensembl's chain file
- Successfully mapped coordinates → `ensembl_hg38entrez_protChange.bed` 
- Still unmapped coordinates → `ensembl_unmapped_entrez_protChange.bed`

### 3. Chromosome Name Standardization
```bash
# Add 'chr' prefix to Ensembl results (chr1, chr2, etc.)
sed 's/^\([a-zA-Z0-9]*\)/chr\1/' ensembl_hg38entrez_protChange.bed
```

### 4. Data Consolidation
- Combines all GRCh38 coordinate files:
  - Original GRCh38 records (`cbio_hg38.bed`)
  - UCSC liftover results (`ucsc_hg38entrez_protChange.bed`)
  - Ensembl liftover results (`ensembl_hg38entrez_protChange.bed`)
- Removes duplicate entries
- Adds standardized column headers

## Output Files

### Intermediate Files
- `cbio_hg38.bed` - Records already in GRCh38 format
- `hg19entrez_protChange.bed` - Records requiring coordinate conversion
- `ucsc_hg38entrez_protChange.bed` - Successfully converted via UCSC chain
- `ucsc_unmapped_entrez_protChange.bed` - Failed UCSC conversions
- `ensembl_hg38entrez_protChange.bed` - Successfully converted via Ensembl chain
- `ensembl_unmapped_entrez_protChange.bed` - Failed conversions from both methods

### Final Output
- `hg38_combined.bed` - Consolidated file with all successfully mapped GRCh38 coordinates
  - **Headers**: chr_id, start_pos, end_pos, entrez_gene_id, prot_change
  - **Format**: Tab-separated values
  - **Content**: All unique records with GRCh38 coordinates

## Error Handling
The script includes comprehensive error checking:
- Validates each AWK operation for data extraction
- Confirms successful liftOver executions
- Exits with error code 1 if any critical step fails
- Provides descriptive error messages for troubleshooting

## Key Features
- **Multi-strategy approach**: Uses two different chain files to maximize conversion success
- **Preserves existing data**: Records already in GRCh38 are maintained without modification
- **Standardized output**: Ensures consistent chromosome naming conventions
- **Duplicate removal**: Eliminates redundant entries in final output
- **Robust error handling**: Comprehensive validation of each processing step

## Dependencies
- bash shell environment
- AWK text processing utility
- `liftOver` tool (UCSC Genome Browser utilities)
- `sed` text stream editor
- Standard Unix utilities (cat, sort, mv)

## Working Directory
All operations are performed in: `/data/shared/biomuta/generated/datasets/2024_10_22/liftover`