import json

# Mapping of DOIDs to predefined cancer types
DOID_TO_CANCER_TYPE = {
    "DOID:11054": "Urinary Bladder Cancer",
    "DOID:1612": "Breast Cancer",
    "DOID:9256": "Colorectal",
    "DOID:5041": "Esophageal Cancer",
    "DOID:11934": "Head and Neck Cancer",
    "DOID:263": "Kidney Cancer",
    "DOID:3571": "Liver Cancer",
    "DOID:1324": "Lung Cancer",
    "DOID:10283": "Prostate Cancer",
    "DOID:10534": "Stomach Cancer",
    "DOID:1781": "Thyroid Gland Cancer",
    "DOID:363": "Uterine Cancer",
    "DOID:4362": "Cervical Cancer",
    "DOID:1319": "Brain Cancer",
    "DOID:2531": "Hematologic Cancer",
    "DOID:3953": "Adrenal Gland Cancer",
    "DOID:1793": "Pancreatic Cancer",
    "DOID:2394": "Ovarian Cancer",
    "DOID:4159": "Skin Cancer"
}

# Mapping of predefined cancer types to TCGA abbreviations
CANCER_TYPE_TO_TCGA = {
    "Urinary Bladder Cancer": "BLCA",
    "Breast Cancer": "BRCA",
    "Colorectal": "COAD",
    "Esophageal Cancer": "ESCA",
    "Head and Neck Cancer": "HNSC",
    "Kidney Cancer": "KICH",
    "Liver Cancer": "LIHC",
    "Lung Cancer": "LUAD",
    "Prostate Cancer": "PRAD",
    "Stomach Cancer": "STAD",
    "Thyroid Gland Cancer": "THCA",
    "Uterine Cancer": "UCEC",
    "Cervical Cancer": "CESC",
    "Brain Cancer": "GBM",
    "Hematologic Cancer": "LAML",
    "Adrenal Gland Cancer": "ACC",
    "Pancreatic Cancer": "PAAD",
    "Ovarian Cancer": "OV",
    "Skin Cancer": "SKCM"
}

def classify_studies(study_mapping_json):
    classified_studies = []
    
    for study in study_mapping_json:
        study_id = study["studyId"]
        do_info = study["do_name"].split(" / ")  # Extract DOID
        doid = do_info[0]
        
        cancer_type = DOID_TO_CANCER_TYPE.get(doid)
        if cancer_type:
            tcga_abbreviation = CANCER_TYPE_TO_TCGA.get(cancer_type, "Unknown")
            classified_studies.append({
                "studyId": study_id,
                "doid": doid,
                "cancerType": cancer_type,
                "TCGA_Study_Abbreviation": tcga_abbreviation
            })
    
    return classified_studies

# Load input JSON
def main():
    with open("/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/study_ids_with_do.json", "r") as f:
        study_mapping_json = json.load(f)
    
    classified_data = classify_studies(study_mapping_json)
    
    with open("/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/classified_studies.json", "w") as f:
        json.dump(classified_data, f, indent=4)
    
    print("Classification complete. Output saved to classified_studies.json")

if __name__ == "__main__":
    main()