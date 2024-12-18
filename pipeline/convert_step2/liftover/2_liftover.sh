#!/bin/bash

wd="/data/shared/biomuta/generated/datasets/2024_10_22/liftover"

# Extract rows with GRCh38 and save as tab-separated:
awk '$5 == "GRCh38"' ${wd}/hg19entrez_build_protChange.bed | awk '{OFS="\t"; $5=""; $1=$1; print}' > ${wd}/cbio_hg38.bed
# Check if the extraction and modification were successful
if [ $? -ne 0 ]; then
    echo "Error extracting or modifying rows for GRCh38."
    exit 1
fi

# Extract all other rows (where the 5th column is not GRCh38) and save as tab-separated:
awk '$5 != "GRCh38"' ${wd}/hg19entrez_build_protChange.bed | awk '{OFS="\t"; $5=""; $1=$1; print}' > ${wd}/hg19entrez_protChange.bed
# Check if the extraction and modification were successful
if [ $? -ne 0 ]; then
    echo "Error extracting or modifying non-GRCh38 rows."
    exit 1
fi

# Run liftOver to convert coordinates (first chain)
./liftOver ${wd}/hg19entrez_protChange.bed ucscHg19ToHg38.over.chain ${wd}/ucsc_hg38entrez_protChange.bed ${wd}/ucsc_unmapped_entrez_protChange.bed

# Check if the first liftOver was successful
if [ $? -ne 0 ]; then
    echo "Error running liftOver with ucscHg19ToHg38.chain."
    exit 1
fi

# Run liftOver to convert coordinates (second chain)
./liftOver ${wd}/ucsc_unmapped_entrez_protChange.bed ensembl_GRCh37_to_GRCh38.chain ${wd}/ensembl_hg38entrez_protChange.bed ${wd}/ensembl_unmapped_entrez_protChange.bed

# Check if the second liftOver was successful
if [ $? -ne 0 ]; then
    echo "Error running liftOver with ensembl_GRCh37_to_GRCh38.chain."
    exit 1
fi

# Prepend 'chr' to the 1st column of ensembl_hg38entrez_protChange.bed
sed 's/^\([a-zA-Z0-9]*\)/chr\1/' ${wd}/ensembl_hg38entrez_protChange.bed > ${wd}/temp && mv ${wd}/temp ${wd}/ensembl_hg38entrez_protChange.bed

# Combine all hg38 files
cat ${wd}/cbio_hg38.bed ${wd}/ucsc_hg38entrez_protChange.bed ${wd}/ensembl_hg38entrez_protChange.bed > ${wd}/hg38_combined.bed
# Remove duplicate rows taking into account extra tabs
awk -v OFS='\t' '{$1=$1; print}' ${wd}/hg38_combined.bed | sort -u > ${wd}/temp && mv ${wd}/temp ${wd}/hg38_combined.bed
# Add headers with tab separation
sed -i '1i chr_id\tstart_pos\tend_pos\tentrez_gene_id\tprot_change' ${wd}/hg38_combined.bed

echo "Script completed successfully."


#chromosome pos GRCh37 -> liftover to GRCh38
#GRCh38 positions -> 
#chr start end entrez prot_change build ID
#7   12222 12222 original               1
#7   12345 12345 after liftover         1