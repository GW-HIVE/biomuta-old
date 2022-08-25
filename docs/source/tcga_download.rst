**Download: TCGA**
------------------

Annotated variant files are downloaded from the ISB-CGC Big Query repository.


**Fields** 
^^^^^^^^^^
Field descriptions for Big Query output available in field_names_descriptions.csv

Additional field descriptions available on `GDC docs <https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/>`_


**Studies**
^^^^^^^^^^^
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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For complete documentation, see the `ISB-CGC Read the Docs pages <https://isb-cancer-genomics-cloud.readthedocs.io/en/latest/>`_

**Step 1 - Gain All Access Requiremenets**
""""""""""""""""""""""""""""""""""""""""""""

Contact Dr. Fabian Seidle and ask for access to the ISB-CGC Big Query repository
    - Example: For the run in Spring-Summer 2022 my (Ned's) personal gwu account was added to the project 'isb-cgc-training-001'
    - All users have up to 1 TB of downloads free, for our purposes we are well under this limit so should not need to pay

Gain access to dbGaP data
    - Apply for access to controlled data at `this website <https://dbgap.ncbi.nlm.nih.gov/aa/wga.cgi?page=login>`_
    - You will need to be approved by a PI that already has access to dbGaP controlled data

For further information see the `ISB-CGC documentation on gaining access <https://isb-cancer-genomics-cloud.readthedocs.io/en/latest/sections/Gaining-Access-To-Controlled-Access-Data.html>`_

**Step 2 - Run downloader R script using R Studio**
""""""""""""""""""""""""""""""""""""""""""""""""""""

**TCGA_mutation_download.R**

Run each line one after the other, instead of the whole script together

Running `library(bigrquery)` and calling this library with `bq_project_query()` (later in the script) will open a browser to login with google credentials
    - Use the google account registered with Fabian for a ISB-CGC project and with dbGAP authorization
    - After logging in, a token will be saved so that you can login through R studio instead

This script will download all mutation data for TCGA. 

There were issues in running this script because the downloaded file was so large. 

In this case run the following scripts in the folder `mutation_download_subscripts`:
    - TCGA_mutation_download_part1.R
    - TCGA_mutation_download_part2.R
    - TCGA_mutation_download_part3.R
    - TCGA_mutation_download_part4.R

These scripts will download a set of the TCGA studies, so that the downloaded file size is smaller. 


**Additional Information**
^^^^^^^^^^^^^^^^^^^^^^^^^^

Go to `MyBinder <https://mybinder.org/>`_ 

For 'Github repository name or URL' enter https://github.com/isb-cgc/ISB-CGC-Demos, then click 'Launch'.

The methods in this tutorial were used to generate the R scripts used to download the data.

**get_field_names.R**

Download a list of all field names for the mutation data, many fields are excluded in the mutation downloader script.

**TCGA_clinical_info_download.R**

Download clinical information for all patients included in the mutation file download. 

**get_field_names_clinical_info.R**

Download a list of all field names for the corresponding clinical data, many fields are excluded in the clinical information downloader script. 