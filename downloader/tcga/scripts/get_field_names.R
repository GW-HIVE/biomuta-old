library(bigrquery)

# Establish what project to access 
project <- "isb-cgc-training-001"

# SQL to gather all field names from the mutation dataset
mut_column_sql <- "SELECT 
    column_name 
FROM isb-cgc-bq.TCGA.INFORMATION_SCHEMA.COLUMNS 
WHERE table_name = 'somatic_mutation_hg38_gdc_current'"

# Launch the query and download field names
mut_column_query <- bq_project_query(project, query=mut_column_sql)

# Create a data frame from the downloaded field names
mut_column_df <- as.vector(unlist(bq_table_download(mut_column_query)))

# Create a csv file with the field names
write.csv(mut_column_df, "TCGA_mutation_field_names.csv", row.names = FALSE)
