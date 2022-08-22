BioMuta Pipeline v-5.0 README
===================================

Last Updated May 2022 by Ned Cauley

Description
-----------

The Biomuta pipeline gathers mutation data from various sources and combines them into a single dataset under common field structure. 

The sources included in BioMuta are: 

- **The Cancer Genome Atlas (TCGA)**

BioMuta gathers mutation data for the following cancers:

- **Urinary Bladder Cancer (DOID:11054)**
- **Breast Cancer (DOID:1612)**
- **Colorectal (DOID:9256)**
- **Esophageal Cancer (DOID:5041)**
- **Head and Neck Cancer (DOID:11934)**
- **Kidney Cancer (DOID:263)**
- **Liver Cancer (DOID:3571)**
- **Lung Cancer (DOID:1324)**
- **Prostate Cancer (DOID:10283)**
- **Stomach Cancer (DOID:10534)**
- **Thyroid Gland Cancer (DOID:1781)**
- **Uterine Cancer (DOID:363)**
- **Cervical Cancer (DOID:4362)**
- **Brain Cancer (DOID:1319)**
- **Hematologic Cancer (DOID:2531)**
- **Adrenal Gland Cancer (DOID:3953)**
- **Pancreatic Cancer (DOID:1793)**
- **Ovarian Cancer (DOID:2394)**
- **Skin Cancer (DOID:4159)**

Running the Pipeline
--------------------

To run the BioMuta pipeine, download the scripts from the HIVE Lab github repo: 
`GW HIVE BioMuta Repository <https://github.com/GW-HIVE/biomuta>`_

Pipeline Overview
-----------------

**Step 1: Download**

In the downloader step, mutation lists will be downloaded from each source. Refer to each individual source below for downloading instructions.

Index for download step:

.. toctree::
   download

**Step 2: Convert**

In the convert step, all resources are formatted to the Biomuta standard for both data and field structure.

Index for convert step:

.. toctree::
    convert

**Step 3: Combine**

In the combined step, all resources are combined into a master dataset.

Index for combined step:

.. toctree::
    combine
   
