**Convert: COSMIC**
====================

**Scripts**
-----------

map_cosmic_tsv.py

**Procedure**
-------------

*Run map_cosmic_tsv.py*
^^^^^^^^^^^^^^^^^^^^^^^

*Summary*
"""""""""

The python script map_cosmic_tsv.py will take the output of the TCGA download step and:
    - Map the data to: 
          - uniprot accessions
          - doid parent terms
    - Rename fields
    - Reformat fields
        - amino acid change and position
        - chromosome id
        - genomic location
        - nucleotide change

*Script Specifications*
"""""""""""""""""""""""

**Running map_cosmic_tsv.py**
#############################

The script must be called from the command line and takes specific command line arguments

Input
#####
    * -c : A path to the cosmic tsv mutation file
    * -m : A path to the folder containing mapping files
    * -d : The name of the doid to cosmic cancer type mapping file
    * -e : The name of the enst to uniprot accession mapping file
    * -o : A path to the the folder to export the final mapped mutations


Output
######
    * A mutation file with COSMIC mutations mapped to doid terms and uniprot accessions

Usage
#####
    * map_cosmic_tsv -h

    *Gives a description of the neccessary commands

    * python map_cosmic_tsv.py -c <path/cosmic_file_name.tsv> -m <path/mapping_folder> -d <doid_mapping_file_name> -e <enst_mapping_file_name> -o <path/output_folder>

    *Runs the script with the given input tsv and outputs a tsv with Biomuta formatting.

'''
*Additional Notes*
------------------

**Mapping**
###########

All the mapping files are available in the scripts repository in the folder pipeline/convert_step2/mapping

The mapping files used for converting the COSMIC tsv are:

**DOID:** cosmic_doid_mapping.csv

COSMIC tissue site terms were mapped to DOID parent terms using the following table:

+---------------------------------------------+--------------------------------------------+
| Primary Site                                | Top_Level_Organ_system                     |
+=============================================+============================================+
| NS                                          | NA                                         |
+---------------------------------------------+--------------------------------------------+
| adrenal_gland                               | DOID:3953 / adrenal gland cancer           |
+---------------------------------------------+--------------------------------------------+
| autonomic_ganglia                           | NA                                         |
+---------------------------------------------+--------------------------------------------+
| biliary_tract                               | DOID:4606 / bile duct cancer               |
+---------------------------------------------+--------------------------------------------+
| bone                                        | DOID:184 / bone cancer                     |
+---------------------------------------------+--------------------------------------------+
| breast                                      | DOID:1612 / breast cancer                  |
+---------------------------------------------+--------------------------------------------+
| central_nervous_system                      | DOID:1319 / brain cancer                   |
+---------------------------------------------+--------------------------------------------+
| cervix                                      | DOID:4362 / cervical cancer                |
+---------------------------------------------+--------------------------------------------+
| endometrium                                 | DOID:363 / uterine cancer                  |
+---------------------------------------------+--------------------------------------------+
| eye                                         | DOID:2174 / ocular cancer                  |
+---------------------------------------------+--------------------------------------------+
| fallopian_tube                              | DOID:1964 / fallopian tube cancer          |
+---------------------------------------------+--------------------------------------------+
| female_genital_tract_(site_indeterminate)   |                                            |
+---------------------------------------------+--------------------------------------------+
| female_genitourinary_system                 | NA                                         |
+---------------------------------------------+--------------------------------------------+
| gastrointestinal_tract_(site_indeterminate) | DOID:3119 / gastrointestinal system cancer |
+---------------------------------------------+--------------------------------------------+
| genital_tract                               | NA                                         |
+---------------------------------------------+--------------------------------------------+
| haematopoietic_and_lymphoid_tissue          | DOID:2531 / hematologic cancer             |
+---------------------------------------------+--------------------------------------------+
| kidney                                      | DOID:263 / kidney cancer                   |
+---------------------------------------------+--------------------------------------------+
| large_intestine                             | DOID:9256 / colorectal cancer              |
+---------------------------------------------+--------------------------------------------+
| liver                                       | DOID:3571 / liver cancer                   |
+---------------------------------------------+--------------------------------------------+
| lung                                        | DOID:1324 / lung cancer                    |
+---------------------------------------------+--------------------------------------------+
| mediastinum                                 | NA                                         |
+---------------------------------------------+--------------------------------------------+
| meninges                                    | DOID:3565 / meningioma                     |
+---------------------------------------------+--------------------------------------------+
| oesophagus                                  | DOID:5041 / esophageal cancer              |
+---------------------------------------------+--------------------------------------------+
| ovary                                       | DOID:2394 / ovarian cancer                 |
+---------------------------------------------+--------------------------------------------+
| pancreas                                    | DOID:1793 / pancreatic cancer              |
+---------------------------------------------+--------------------------------------------+
| paratesticular_tissues                      | NA                                         |
+---------------------------------------------+--------------------------------------------+
| parathyroid                                 | DOID:1540 / parathyroid carcinoma          |
+---------------------------------------------+--------------------------------------------+
| penis                                       | DOID:11615 / penile cancer                 |
+---------------------------------------------+--------------------------------------------+
| pericardium                                 | NA                                         |
+---------------------------------------------+--------------------------------------------+
| perineum                                    | DOID:4045 / muscle cancer                  |
+---------------------------------------------+--------------------------------------------+
| peritoneum                                  | DOID:1725 / peritoneum cancer              |
+---------------------------------------------+--------------------------------------------+
| pituitary                                   | DOID:1785 / pituitary cancer               |
+---------------------------------------------+--------------------------------------------+
| placenta                                    | DOID:2021 / placenta cancer                |
+---------------------------------------------+--------------------------------------------+
| pleura                                      | DOID:5158 / pleural cancer                 |
+---------------------------------------------+--------------------------------------------+
| prostate                                    | DOID:10283 / prostate cancer               |
+---------------------------------------------+--------------------------------------------+
| retroperitoneum                             | DOID:5875 / retroperitoneal cancer         |
+---------------------------------------------+--------------------------------------------+
| salivary_gland                              | DOID:8618 / oral cavity cancer             |
+---------------------------------------------+--------------------------------------------+
| skin                                        | DOID:4159 / skin cancer                    |
+---------------------------------------------+--------------------------------------------+
| small_intestine                             | DOID:9253 / gastrointestinal stromal tumor |
+---------------------------------------------+--------------------------------------------+
| soft_tissue                                 | NA                                         |
+---------------------------------------------+--------------------------------------------+
| stomach                                     | DOID:10534 / stomach cancer                |
+---------------------------------------------+--------------------------------------------+
| testis                                      | DOID:2998 / testicular cancer              |
+---------------------------------------------+--------------------------------------------+
| thymus                                      | DOID:3277 / thymus cancer                  |
+---------------------------------------------+--------------------------------------------+
| thyroid                                     | DOID:1781 / thyroid gland cancer           |
+---------------------------------------------+--------------------------------------------+
| upper_aerodigestive_tract                   | DOID:8618 / oral cavity cancer             |
+---------------------------------------------+--------------------------------------------+
| urinary_tract                               | DOID:11054 / urinary bladder cancer        |
+---------------------------------------------+--------------------------------------------+
| uterine_adnexa                              | NA                                         |
+---------------------------------------------+--------------------------------------------+
| vagina                                      | DOID:119 / vaginal cancer                  |
+---------------------------------------------+--------------------------------------------+
| vulva                                       | DOID:1245 / vulva cancer                   |
+---------------------------------------------+--------------------------------------------+

**Uniprot Accession:** human_protein_transcriptlocus.csv

Transcript ID (starts with ENST) was mapped to uniprot isoform accession

Mapping was NOT performed to uniprot canonical accession as this resulted in an issue with the final dataset in which a mutation for the same canonical accession would be listed with different amino acid changes