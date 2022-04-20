**************
BioMuta README
**************

Updated April 2022 by Ned Cauley

Overview
########

BioMuta is a pipeline for generating a single-nucleotide variation (SNV) and disease association dataset where variations are mapped to genomes and RefSeq nucleotide entries, and unified through UniProtKB/Swiss-Prot positional coordinates.

BioMuta gathers SNV data for the following cancers:

* Urinary Bladder Cancer (DOID:11054)
* Breast Cancer (DOID:1612)
* Colorectal (DOID:9256)
* Esophageal Cancer (DOID:5041)
* Head and Neck Cancer (DOID:11934)
* Kidney Cancer (DOID:263)
* Liver Cancer (DOID:3571)
* Lung Cancer (DOID:1324)
* Prostate Cancer (DOID:10283)
* Stomach Cancer (DOID:10534)
* Thyroid Gland Cancer (DOID:1781)
* Uterine Cancer (DOID:363)
* Cervical Cancer (DOID:4362)
* Brain Cancer (DOID:1319)
* Hematologic Cancer (DOID:2531)
* Adrenal Gland Cancer (DOID:3953)
* Pancreatic Cancer (DOID:1793)
* Ovarian Cancer (DOID:2394)
* Skin Cancer (DOID:4159)

Running the Pipeline
####################

To run the BioMuta pipeline, you will need to download the required scripts and mapping files from the `BioMuta GW HIVE lab Github repository <https://github.com/GW-HIVE/biomuta>`_

Pipeline Overview
*****************

Step 1: Downloader
******************

Resource: TCGA
--------------

Annotated variant files are downloaded from the ISB-CGC Big Query repository.


**Fields** 
Field descriptions for Big Query output available in field_names_descriptions.csv

Additional field descriptions available on `GDC docs <https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/>`


**Studies**

+---------------+------------------------------------------------------------------+
| TCGA Study ID | TCGA Study Name                                                  |
+===============+==================================================================+
| ACC           | Adrenocortical carcinoma                                         |
+---------------+------------------------------------------------------------------+
| BLCA          | Bladder Urothelial Carcinoma                                     |
+---------------+------------------------------------------------------------------+
| BRCA          | Breast invasive carcinoma                                        |
+---------------+------------------------------------------------------------------+
| CESC          | Cervical squamous cell carcinoma and endocervical adenocarcinoma |
+---------------+------------------------------------------------------------------+
| CHOL          | Cholangiocarcinoma                                               |
+---------------+------------------------------------------------------------------+
| COAD          | Colon adenocarcinoma                                             |
+---------------+------------------------------------------------------------------+
| DLBC          | Lymphoid Neoplasm Diffuse Large B-cell Lymphoma                  |
+---------------+------------------------------------------------------------------+
| ESCA          | Esophageal carcinoma                                             |
+---------------+------------------------------------------------------------------+
| GBM           | Glioblastoma multiforme                                          |
+---------------+------------------------------------------------------------------+
| HNSC          | Head and Neck squamous cell carcinoma                            |
+---------------+------------------------------------------------------------------+
| KICH          | Kidney Chromophobe                                               |
+---------------+------------------------------------------------------------------+
| KIRC          | Kidney renal clear cell carcinoma                                |
+---------------+------------------------------------------------------------------+
| KIRP          | Kidney renal papillary cell carcinoma                            |
+---------------+------------------------------------------------------------------+
| LAML          | Acute Myeloid Leukemia                                           |
+---------------+------------------------------------------------------------------+
| LGG           | Brain Lower Grade Glioma                                         |
+---------------+------------------------------------------------------------------+
| LIHC          | Liver hepatocellular carcinoma                                   |
+---------------+------------------------------------------------------------------+
| LUAD          | Lung adenocarcinoma                                              |
+---------------+------------------------------------------------------------------+
| LUSC          | Lung squamous cell carcinoma                                     |
+---------------+------------------------------------------------------------------+
| MESO          | Mesothelioma                                                     |
+---------------+------------------------------------------------------------------+
| OV            | Ovarian serous cystadenocarcinoma                                |
+---------------+------------------------------------------------------------------+
| PAAD          | Pancreatic adenocarcinoma                                        |
+---------------+------------------------------------------------------------------+
| PCPG          | Pheochromocytoma and Paraganglioma                               |
+---------------+------------------------------------------------------------------+
| PRAD          | Prostate adenocarcinoma                                          |
+---------------+------------------------------------------------------------------+
| READ          | Rectum adenocarcinoma                                            |
+---------------+------------------------------------------------------------------+
| SARC          | Sarcoma                                                          |
+---------------+------------------------------------------------------------------+
| SKCM          | Skin Cutaneous Melanoma                                          |
+---------------+------------------------------------------------------------------+
| STAD          | Stomach adenocarcinoma                                           |
+---------------+------------------------------------------------------------------+
| TGCT          | Testicular Germ Cell Tumors                                      |
+---------------+------------------------------------------------------------------+
| THCA          | Thyroid carcinoma                                                |
+---------------+------------------------------------------------------------------+
| THYM          | Thymoma                                                          |
+---------------+------------------------------------------------------------------+
| UCEC          | Uterine Corpus Endometrial Carcinoma                             |
+---------------+------------------------------------------------------------------+
| UCS           | Uterine Carcinosarcoma                                           |
+---------------+------------------------------------------------------------------+
| UVM           | Uveal Melanoma                                                   |
+---------------+------------------------------------------------------------------+

**Downloading through Big Query**
*********************************

*Gain access to the ISB-GCG Big Query repository*

Contact Dr. Fabian Seidle and ask for access. 
    - For the run in Spring 2022 my (Ned) personal gwu account was added to the project 'isb-cgc-training-001'
    - All users have up to 1 TB of downloads free, for our ourporposes we are well under this limit so should not need to pay



*Step 2 - Review the tutorial offered by the ISB-CGC*

Go to `MyBinder <https://mybinder.org/>` 

For 'Github repository name or URL' enter https://github.com/isb-cgc/ISB-CGC-Demos, then click 'Launch'.

The methods in this tutorial were used to generate the R scripts used to download the data.










Version 4.0 Information (deprecated)
####################################

Previous versions' script for automatic download of source files is deprecated due to several resources changing access methods.

Resource: TCGA
--------------

VCF files are generated through the GDC through matching tumor and normal samples as detailed `here <https://docs.gdc.cancer.gov/Data/File_Formats/VCF_Format/>`_

**TCGA Download**

Go to the `GDC data portal repository <https://portal.gdc.cancer.gov/repository>`_

Go to `Advanced Search` and enter the following query:

``cases.project.program.name in ["TCGA"] and cases.project.project_id in ["TCGA-ACC","TCGA-BLCA","TCGA-BRCA","TCGA-CESC","TCGA-COAD","TCGA-DLBC","TCGA-ESCA","TCGA-GBM","TCGA-HNSC","TCGA-KICH","TCGA-KIRC","TCGA-KIRP","TCGA-LAML","TCGA-LGG","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-OV","TCGA-PAAD","TCGA-PCPG","TCGA-PRAD","TCGA-READ","TCGA-SKCM","TCGA-STAD","TCGA-THCA","TCGA-UCEC","TCGA-UCS"] and files.data_format in ["vcf"] and files.data_type in ["Raw Simple Somatic Mutation"]``

This query will yield 5 times as many files as there are cases. This is because each case was run through 5 separate variant calling pipelines.

The variant calling pipelines are detailed `here <https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/#somatic-variant-calling-workflow>`_

The most recent variant calling tool incorporated is `MuSe <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1029-6>`_

After identifying workflow type, the following still have more files then cases: 
* TCGA-BLCA
* TCGA-BRCA
* TCGA-CESC
* and others...

Why is this the case?

Investigation of TCGA-ESCA:

Found that the repeated sample is TCGA-V5-A7RC
    * On GDC data portal, TCGA-V5-A7RC contained one sample for primary tumor (similar to other 183 cases) and one sample for metastatic (unique)

We can obtain only unique samples by using the search: 

`cases.project.program.name in ["TCGA"] and cases.project.project_id in ["TCGA-ESCA"] and files.analysis.workflow_type in ["MuSE"] and files.data_format in ["vcf"] and cases.samples.sample_type in ["primary tumor"]`

Do we eliminate this metastatic sample? Or do we include samples from the same case? Will it matter to have multiple samples form the same case downstream?

For TCGA-BRCA: 

Updating the search to include primary tumor only does not eliminate all duplicates. 
    * Still found 7 metastatic samples

For specific case, three duplicates but no differences between samples identified:
    * `https://portal.gdc.cancer.gov/repository?facetTab=cases&filters=%7B%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.case_id%22%2C%22value%22%3A%5B%22f130f376-5801-40f9-975d-a7e2f7b5670d%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.program.name%22%2C%22value%22%3A%5B%22TCGA%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22TCGA-BRCA%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.analysis.workflow_type%22%2C%22value%22%3A%5B%22MuSE%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_format%22%2C%22value%22%3A%5B%22vcf%22%5D%7D%7D%5D%2C%22op%22%3A%22and%22%7D`

Selected 

Repeat this query for all desired TCGA studies. The following studies were gathered for v5.0:

* BLCA
* BRCA
* COAD
* ESCA
* HNSC
* KICH
* KIRC
* KIRP
* LIHC
* LUAD
* LUSC
* PRAD
* READ
* STAD
* THCA
* UCEC

* Pancreatic Cancer
    * PAAD
* Ovarian Cancer
    * OV
* Skin Cancer
    * SKCM
* Cervical Cancer
    * CESC
* Brain Cancer
    * LGG
    * GBM
* Hematologic Cancer
    * LAML
    * DLBC
* Adrenal Gland Cancer
    * ACC
    * PCPG
* Uterine Cancer
    * UCS

TCGA studies from v-4.0 removed in v-5.0

* CHOL - Cholangiocarcinoma
    * Bile Duct Cancer
* MESO - Mesothelioma
    * Malignant Mesothelioma
* SARC - Sarcoma
    *  Sarcoma Cell
* TGCT - Testicular Germ Cell Tumors
    * Testicular Cancer
* THYM - Thymoma
    * Thymus Cancer
* UVM - Uveal Melanoma
    * Ocular Cancer


**Step 1: Downloader**

In the downloader step, SNV data is downloaded directly from sources using either a manual download from the source's web page or a call to an API.

**Step 2: Annotator**

Step 1: Downloader
******************
Previous versions' script for automatic download of source files is deprecated due to several resources changing access methods.

Resource: TCGA
--------------

VCF files are generated through the GDC through matching tumor and normal samples as detailed `here <https://docs.gdc.cancer.gov/Data/File_Formats/VCF_Format/>`_

**TCGA Download**

Go to the `GDC data portal repository <https://portal.gdc.cancer.gov/repository>`_

Go to `Advanced Search` and enter the following query:

``cases.project.program.name in ["TCGA"] and cases.project.project_id in ["TCGA-ACC","TCGA-BLCA","TCGA-BRCA","TCGA-CESC","TCGA-COAD","TCGA-DLBC","TCGA-ESCA","TCGA-GBM","TCGA-HNSC","TCGA-KICH","TCGA-KIRC","TCGA-KIRP","TCGA-LAML","TCGA-LGG","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-OV","TCGA-PAAD","TCGA-PCPG","TCGA-PRAD","TCGA-READ","TCGA-SKCM","TCGA-STAD","TCGA-THCA","TCGA-UCEC","TCGA-UCS"] and files.data_format in ["vcf"] and files.data_type in ["Raw Simple Somatic Mutation"]``

This query will yield 5 times as many files as there are cases. This is because each case was run through 5 separate variant calling pipelines.

The variant calling pipelines are detailed `here <https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/#somatic-variant-calling-workflow>`_

The most recent variant calling tool incorporated is `MuSe <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1029-6>`_

After identifying workflow type, the following still have more files then cases: 
* TCGA-BLCA
* TCGA-BRCA
* TCGA-CESC
* and others...

Why is this the case?

Investigation of TCGA-ESCA:

Found that the repeated sample is TCGA-V5-A7RC
    * On GDC data portal, TCGA-V5-A7RC contained one sample for primary tumor (similar to other 183 cases) and one sample for metastatic (unique)

We can obtain only unique samples by using the search: 

`cases.project.program.name in ["TCGA"] and cases.project.project_id in ["TCGA-ESCA"] and files.analysis.workflow_type in ["MuSE"] and files.data_format in ["vcf"] and cases.samples.sample_type in ["primary tumor"]`

Do we eliminate this metastatic sample? Or do we include samples from the same case? Will it matter to have multiple samples form the same case downstream?

For TCGA-BRCA: 

Updating the search to include primary tumor only does not eliminate all duplicates. 
    * Still found 7 metastatic samples

For specific case, three duplicates but no differences between samples identified:
    * `https://portal.gdc.cancer.gov/repository?facetTab=cases&filters=%7B%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.case_id%22%2C%22value%22%3A%5B%22f130f376-5801-40f9-975d-a7e2f7b5670d%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.program.name%22%2C%22value%22%3A%5B%22TCGA%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22TCGA-BRCA%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.analysis.workflow_type%22%2C%22value%22%3A%5B%22MuSE%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_format%22%2C%22value%22%3A%5B%22vcf%22%5D%7D%7D%5D%2C%22op%22%3A%22and%22%7D`

Selected 

Repeat this query for all desired TCGA studies. The following studies were gathered for v5.0:

* BLCA
* BRCA
* COAD
* ESCA
* HNSC
* KICH
* KIRC
* KIRP
* LIHC
* LUAD
* LUSC
* PRAD
* READ
* STAD
* THCA
* UCEC

* Pancreatic Cancer
    * PAAD
* Ovarian Cancer
    * OV
* Skin Cancer
    * SKCM
* Cervical Cancer
    * CESC
* Brain Cancer
    * LGG
    * GBM
* Hematologic Cancer
    * LAML
    * DLBC
* Adrenal Gland Cancer
    * ACC
    * PCPG
* Uterine Cancer
    * UCS

TCGA studies from v-4.0 removed in v-5.0

* CHOL - Cholangiocarcinoma
    * Bile Duct Cancer
* MESO - Mesothelioma
    * Malignant Mesothelioma
* SARC - Sarcoma
    *  Sarcoma Cell
* TGCT - Testicular Germ Cell Tumors
    * Testicular Cancer
* THYM - Thymoma
    * Thymus Cancer
* UVM - Uveal Melanoma
    * Ocular Cancer

Resource: COSMIC
----------------

`COSMIC <https://cancer.sanger.ac.uk/cosmic/download>`_

What are the disease included in COSMIC?

**COSMIC Download**

To download the COSMIC mutation set, run the script ``download_cosmic.py``

    Error: downloading an empty file for some reason...

Resource: ClinVar
-----------------

`ClinVar <https://www.ncbi.nlm.nih.gov/clinvar/>`_

**ClinVar download**

Use the FTP site for access to the latetest clinVar release: 

``wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz``

For more information see the `README <https://ftp.ncbi.nlm.nih.gov/pub/clinvar/README.txt>`_

Resource: ICGC
--------------

`ICGC <https://dcc.icgc.org/>`_

**ICGC Download**

Manual download at the `Data Portal <https://dcc.icgc.org/releases/release_28/Summary>`

For v-5.0: Downloaded file ``simple_somatic_mutation.aggregated.vcf.gz``

Resource: Intogen
-----------------

`IntOGen <https://www.intogen.org/search>`_

**IntOGen Download**

Manual Download from the `downloads page <https://www.intogen.org/download>`_

For v-5.0: Downloaded file ``intogen_driver_mutations_catalog-2016.5.gz``


Big Query Investigation
-----------------------

Datasets available through Big Query from ISB

List of datasets: https://isb-cancer-genomics-cloud.readthedocs.io/en/latest/sections/Hosted-Data.html

