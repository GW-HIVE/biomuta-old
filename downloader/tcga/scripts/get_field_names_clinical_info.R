library(bigrquery)

clinical_column_sql <- "SELECT 
    column_name 
FROM isb-cgc-bq.TCGA.INFORMATION_SCHEMA.COLUMNS 
WHERE table_name = 'clinical_gdc_current'"

clinical_column_query <- bq_project_query

clinical_case_df <- as.vector(unlist(bq_table_download(clinical_column_query)))

write.csv(mut_df, "clinical_information_fields.csv", row.names = FALSE)