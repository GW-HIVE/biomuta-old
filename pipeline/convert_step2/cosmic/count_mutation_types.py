import logging
import pandas as pd

logging.basicConfig(filename="cosmic.log",
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    filemode='w',
                    level=logging.INFO)
logging.info("Logging started ------------------------")

# Load the data into a Pandas DataFrame
file_path = "uniprot_export.tsv"
df = pd.read_csv(file_path, sep="\t", dtype={"CHROMOSOME": str})
logging.info("Successfully read the TSV file into a pandas dataframe.")

# Initialize counters for mutation types
counts = {
    "SNP": 0,
    "Insertion": 0,
    "Deletion": 0,
    "Frameshift": 0,
    "Indel": 0,
    "Inversion": 0,
    "Others": 0
}

# Define a function to classify mutations
def classify_mutation(row):
    # Single Nucleotide Polymorphism (SNP)
    if not pd.isna(row['GENOMIC_WT_ALLELE_SEQ']) and not pd.isna(row['GENOMIC_MUT_ALLELE_SEQ']):
        if len(row['GENOMIC_WT_ALLELE_SEQ']) == len(row['GENOMIC_MUT_ALLELE_SEQ']) == 1:
            return "SNP"
    # Insertion
    elif pd.isna(row['GENOMIC_WT_ALLELE_SEQ']):
        return "Insertion"
    # Deletion
    elif pd.isna(row['GENOMIC_MUT_ALLELE_SEQ']):
        return "Deletion"
    
    # Frameshift and indels (based on CDS_MUT_SYNTAX or AA_MUT_SYNTAX)
    if "fs" in str(row['AA_MUT_SYNTAX']) or "fs" in str(row['CDS_MUT_SYNTAX']):
        return "Frameshift"
    if "delins" in str(row['AA_MUT_SYNTAX']) or "delins" in str(row['CDS_MUT_SYNTAX']):
        return "Indel"
    if "inv" in str(row['CDS_MUT_SYNTAX']):
        return "Inversion"
    
    # Catch-all for others
    return "Others"

# Apply the classification function to the DataFrame
logging.info("Classifying mutations...")
df['Mutation_Type'] = df.apply(classify_mutation, axis=1)

# Count mutation types
logging.info("Counting mutation types...")
counts.update(df['Mutation_Type'].value_counts().to_dict())

# Print the counts
for mutation_type, count in counts.items():
    logging.info(f"{mutation_type}: {count}")

# Filter rows classified as "SNP" and write to a file
logging.info("Filtering SNPs...")
snps = df[df['Mutation_Type'] == "SNP"]
cosmic_snps = '/data/shared/repos/biomuta-old/downloads/cosmic/cosmic_snps.tsv'
columns_to_omit = ['STRAND', 'COSV_ID', 'CDS_MUT_SYNTAX', 'COSMIC_LEGACY_ID', 'Mutation_Type']
with open(cosmic_snps, 'w') as f:
    snps.drop(columns=columns_to_omit).to_csv(f, sep="\t", index=False)