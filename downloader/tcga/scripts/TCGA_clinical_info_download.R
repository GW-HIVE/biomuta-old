library(bigrquery)

clinical_case_sql <- "SELECT
  submitter_id, 
  primary_site, 
  disease_type, 
  proj__name, 
  proj__project_id, 
  demo__demographic_id, 
  demo__gender, 
  demo__race, 
  demo__ethnicity, 
  demo__vital_status, 
  demo__days_to_birth, 
  demo__year_of_birth, 
  demo__age_at_index, 
  demo__year_of_death, 
  demo__days_to_death, 
  demo__state, 
  demo__created_datetime, 
  demo__updated_datetime, 
  diag__diagnosis_id, 
  diag__ajcc_clinical_n, 
  diag__masaoka_stage, 
  diag__ajcc_clinical_m, 
  diag__primary_diagnosis, 
  diag__primary_gleason_grade, 
  diag__year_of_diagnosis, 
  diag__figo_stage, 
  diag__progression_or_recurrence, 
  diag__ajcc_pathologic_m, 
  diag__site_of_resection_or_biopsy, 
  diag__ajcc_staging_system_edition, 
  diag__icd_10_code, 
  diag__age_at_diagnosis, 
  diag__ajcc_clinical_t, 
  diag__days_to_last_follow_up, 
  diag__ajcc_pathologic_stage, 
  diag__tumor_grade, 
  diag__last_known_disease_status
FROM `isb-cgc-bq.TCGA.clinical_gdc_current`
"

clinical_case_query <- bq_project_query(project, query=clinical_case_sql)

clinical_case_df <- as.data.frame(bq_table_download(clinical_case_query))

write.csv(mut_df, "clinical_information.csv", row.names = FALSE)