library(bigrquery)

# Establish what project to access 
project <- "isb-cgc-training-001"

# SQL to gather mutation data with listed fields
mut_sql <- "SELECT project_short_name, 
    case_barcode, 
    case_id, 
    Hugo_Symbol, 
    Gene, 
    Transcript_ID, 
    ENSP, 
    Variant_Classification, 
    Variant_Type, 
    Chromosome, 
    Start_Position, 
    End_Position, 
    Reference_Allele, 
    Tumor_Seq_Allele1, 
    Tumor_Seq_Allele2, 
    Allele, 
    SIFT, 
    Polyphen, 
    GMAF, 
    EXAC_AF, 
    Amino_acids, 
    CDS_position, 
    Protein_position, 
    dbSNP_RS
FROM `isb-cgc-bq.TCGA.somatic_mutation_hg38_gdc_current`
WHERE Variant_Type = 'SNP' 
    AND Variant_Classification IN ('Missense_Mutation', 
        'Nonsense_Mutation', 
        'Nonstop_Mutation') 
  AND project_short_name IN ('TCGA-GBM', 
        'TCGA-HNSC', 
        'TCGA-KICH', 
        'TCGA-KIRP', 
        'TCGA-KIRC', 
        'TCGA-LAML', 
        'TCGA-LGG', 
        'TCGA-LIHC')"

# Launch the query and download data
mut_query <- bq_project_query(project, query=mut_sql)

# Create a data frame from the downloaded data
mut_df <- bq_table_download(mut_query)

# Create a csv file with the field names
write.csv(mut_df, "TCGA_SNP_somatic_mutation_hg38_part2.csv", row.names = FALSE)