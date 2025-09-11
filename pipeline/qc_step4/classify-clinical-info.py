import csv
import json
import pandas as pd

# Load classified mapping
def load_mapping():
    with open("/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/classified_studies.json", "r") as f:
        classified_studies = json.load(f)
    return {study["studyId"]: f"TCGA-{study['TCGA_Study_Abbreviation']}" for study in classified_studies}

def update_csv(mapping):
    # Load clinical information CSV
    df = pd.read_csv("/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/clinical-information-sample.csv")
    
    # Replace project IDs with TCGA abbreviations
    df["proj__project_id"] = df["proj__project_id"].map(mapping).fillna("Unknown")
    
    # Save updated CSV
    df.to_csv("/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/updated-clinical-information.csv", index=False, quoting=csv.QUOTE_ALL)
    print("CSV updated and saved as updated-clinical-information.csv")

def main():
    mapping = load_mapping()
    update_csv(mapping)

if __name__ == "__main__":
    main()
