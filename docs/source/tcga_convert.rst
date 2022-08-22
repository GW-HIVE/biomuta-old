**Convert: TCGA**
=================

**Scripts**
-----------

process_tcga_download.py

**Procedure**
-------------

*Run process_tcga_download.py*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*Summary*
"""""""""

The python script process_tcga_download.py will take the output of the TCGA download step and:
    - Map the data to: 
          - uniprot accession
          - doid parent terms
    - Rename fields
    - Reformat fields
        - amino acid change and position
        - chromosome id
    - Filter out unnecessary fields 

*Script Specifications*
"""""""""""""""""""""""

**Running process_tcga_download.py**
####################################

The script must be called from the command line and takes specific command line arguments

Input
#####
    * -i : A path to the input csv to reformat
    * -m : A path to the folder containing mapping files
    * -d : A path to the tcga study to doid mapping file
    * -e : A path to the ENSP to uniprot mapping file
    * -o : A path to the output folder

Output
######
    * A data report comparing new AA sites to old AA sites for Biomuta

Usage
#####
    * python process_tcga_download.py -h

    *Gives a description of the neccessary commands

    * python process_tcga_download.py -i <path/input_file.vcf> -m <path/> -d <doid_mapping.csv -e <ensp_mapping.csv> -o <path/>

    *Runs the script with the given input tcga csv and output a formatted csv

*Additional Notes*
------------------

**Mapping**
###########

All the mapping files are available in the scripts repository in the folder pipeline/convert_step2/mapping

The mapping files used for converting TCGA are:

**DOID:** tcga_doid_mapping.csv

TCGA Projects were mapped to DOID parent terms using the following table:

+------------+------------------------+--------------+
| DO_slim_id | DO_slim_name           | TCGA_project |
+============+========================+==============+
| DOID:5041  | esophageal cancer      | TCGA-ESCA    |
+------------+------------------------+--------------+
| DOID:2531  | hematologic cancer     | TCGA-DLBC    |
+------------+------------------------+--------------+
| DOID:9256  | colorectal cancer      | TCGA-READ    |
+------------+------------------------+--------------+
| DOID:1319  | brain cancer           | TCGA-GBM     |
+------------+------------------------+--------------+
| DOID:1319  | brain cancer           | TCGA-LGG     |
+------------+------------------------+--------------+
| DOID:1781  | thyroid cancer         | TCGA-THCA    |
+------------+------------------------+--------------+
| DOID:11054 | urinary bladder cancer | TCGA-BLCA    |
+------------+------------------------+--------------+
| DOID:363   | uterine cancer         | TCGA-UCEC    |
+------------+------------------------+--------------+
| DOID:169   | neuroendocrine tumor   | TCGA-PCPG    |
+------------+------------------------+--------------+
| DOID:4362  | cervical cancer        | TCGA-CESC    |
+------------+------------------------+--------------+
| DOID:363   | uterine cancer         | TCGA-UCS     |
+------------+------------------------+--------------+
| DOID:3277  | thymus cancer          | TCGA-THYM    |
+------------+------------------------+--------------+
| DOID:3571  | liver cancer           | TCGA-LIHC    |
+------------+------------------------+--------------+
| DOID:11934 | head and neck cancer   | TCGA-HNSC    |
+------------+------------------------+--------------+
| DOID:2174  | ocular cancer          | TCGA-UVM     |
+------------+------------------------+--------------+
| DOID:4159  | skin cancer            | TCGA-SKCM    |
+------------+------------------------+--------------+
| DOID:9256  | colorectal cancer      | TCGA-COAD    |
+------------+------------------------+--------------+
| DOID:3953  | adrenal gland cancer   | TCGA-ACC     |
+------------+------------------------+--------------+
| DOID:1793  | pancreatic cancer      | TCGA-PAAD    |
+------------+------------------------+--------------+
| DOID:2994  | germ cell cancer       | TCGA-TGCT    |
+------------+------------------------+--------------+
| DOID:1324  | lung cancer            | TCGA-LUSC    |
+------------+------------------------+--------------+
| DOID:1790  | malignant mesothelioma | TCGA-MESO    |
+------------+------------------------+--------------+
| DOID:2394  | ovarian cancer         | TCGA-OV      |
+------------+------------------------+--------------+
| DOID:1115  | sarcoma                | TCGA-SARC    |
+------------+------------------------+--------------+
| DOID:263   | kidney cancer          | TCGA-KIRP    |
+------------+------------------------+--------------+
| DOID:10534 | stomach cancer         | TCGA-STAD    |
+------------+------------------------+--------------+
| DOID:2531  | hematologic cancer     | TCGA-LAML    |
+------------+------------------------+--------------+
| DOID:10283 | prostate cancer        | TCGA-PRAD    |
+------------+------------------------+--------------+
| DOID:1324  | lung cancer            | TCGA-LUAD    |
+------------+------------------------+--------------+
| DOID:1612  | breast cancer          | TCGA-BRCA    |
+------------+------------------------+--------------+
| DOID:263   | kidney cancer          | TCGA-KIRC    |
+------------+------------------------+--------------+
| DOID:263   | kidney cancer          | TCGA-KICH    |
+------------+------------------------+--------------+

**Uniprot Accession:** human_protein_transcriptlocus.csv

Peptide ID (starts with ENSP) was mapped to uniprot isoform accession

Mapping was NOT performed to uniprot canonical accession as this resulted in an issue with the final dataset in which a mutation for the same canonical accession would be listed with different amino acid changes