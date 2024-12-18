**Convert: ICGC**
=================

**Scripts**
-----------

genomic liftover (mapvcf_copySA.py) > convert_icgc_vcf.py > map_icgc.py

**Procedure**
-------------
*Perform liftover of mutations from GRCh37 to GRCh38 (mapvcf_copySA)*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*Summary*
"""""""""

The most recent data release for CIVIC is aligned to the GRCH37 human reference genome. For this update however, we are using the human reference genome GRCh38.

To convert coordinates between the two reference genomes, we use a 'liftover' tool to remap the genomic coordinates. 

Seun performed the liftover and provided the notes listed below.

**Genomic Liftover Notes**
^^^^^^^^^^^^^^^^^^^^^^^^^^
Step performed and notes provided by Seun Agbaje

**VCF**

A VCF (Variant Call Format) file is a text file used to store gene sequence variations. The files often start with lines of metadata, then headers relating to the variants described. Because the standard for formatting and relaying genomic data is always evolving, there are numerous versions and references for VCF files and the dependencies they use


**Fields**
""""""""""
Common fields for VCF files include:

+---------+----------------------------------------------------+
| Field   | Description                                        |
+=========+====================================================+
| Chrom   | chromosome that the variation is being called on   |
+---------+----------------------------------------------------+
| Pos     | 1 base position of the variant                     |
+---------+----------------------------------------------------+
| ID      | identifier of the variant                          |
+---------+----------------------------------------------------+
| Ref     | reference base at the position of variance         |
+---------+----------------------------------------------------+
| Alt     | alternate alleles at the position                  |
+---------+----------------------------------------------------+
| Qual    | quality score ofthe given alleles                  |
+---------+----------------------------------------------------+
| Filter  | indicates which set of filters failed or passed    |
+---------+----------------------------------------------------+
| Info    | descriptions of the variation                      |
+---------+----------------------------------------------------+
| Format  | (optional) fields describing the sample            |
+---------+----------------------------------------------------+
| Samples | values for each of the samples listed under format |
+---------+----------------------------------------------------+

**Converting with CrossMap**
""""""""""""""""""""""""""""

CrossMap is a program that can convert genome coordinates between different assemblies, such as hg18 (GRCh36) to hg19 (GRCh37). It is made in python and offered as a webtool, by Ensembl in limited capacity or as a local script For full functionality. This gives extra customizability and the option to convert files over 50 mb, it is necessary to run a local edition of CrossMap.

Crossmap Documentation: http://crossmap.sourceforge.net/

*Requirements* 

- Python2 or Python3 installed on a linux server
- Chain file - describes a pairwise alignment between two reference assemblies
- They can be found through UCSC, Ensembl, and other sources
- compressed files are allowed
- hg19ToHg38.over.chain was best tested
- target, input file - file to be converted in format compatible with CrossMap
- CrossMap supports vcf, bam/cram/sam, maf, and other formats torelay genomic data
- compressed files are allowed
- referencefile - fasta format of the wanted genome assembly

*Other files used*

- mapvcf is the script from the package that does the conversion. attachedis the version I used. I believe commenting out lines 100:109 is what allowed it to work
- hg19ToHg38 is the chain file that I used
- this is the command I used to get the assembly file, which is from UCSC "wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz &"
- this is the command I used to unzip the assembly file "gzip -dk hg38.fa.gz"
- the exact command I ran to create the file is  this "python3 ./.local/bin/CrossMap.py vcf /mnt/d/hg19ToHg38.over.chain.gz /mnt/d/icgc_missense_mutations.vcf hg38.fa /mnt/d/icgc_missense_mutations_38_hgfz.vcf"
- of note, there are numerous other assembly and chain files. I tried 3 or 4 of each and the ones linked here were the best. I determined best by both what the script relays and how big the final vcf file were

**Output**

Two output files were generated from the liftover and stored on the OncoMX-tst server at /software/pipeline/integrator/downloads/biomuta/v-5.0/icgc/
- icgc_missense_mutations_38.vcf
  - All mutations with converted coordinates
- icgc_missense_mutations_38_fail.vcf
  - Mutations whose coordinates could not be converted

Only the mutations whose coordinates were successfully converted were carried forward in the pipeline

*Run convert_icgc_vcf.py*
^^^^^^^^^^^^^^^^^^^^^^^^^

*Summary*
"""""""""

The python script convert_icgc_vcf.py will convert the vcf formatted mutation file to a csv file.

With the vcf format, each mutation line in the file can contain multiple annotations and annotation-specific information.

The output csv format will contain only one annotation per line with associated annotation-sepcific information. 

In order to know how the information for the mutation and annotation fields are structured, a schema describing the fields is provided to the script.

**Example Line Tranformation**

*Input VCF lines*

mutation A info | mutation A annotation 1 info | mutation A annotation 2 info

mutation B info | mutation B annotation 1 info | mutation B annotation 2 info | mutation B annotation 3 info

*Output CSV lines*

mutation A info,annotation 1 info

mutation A info,annotation 2 info

mutation B info,annotation 1 info

mutation B info,annotation 2 info

mutation B info,annotation 3 info

*Script Specifications*
"""""""""""""""""""""""

The script must be called from the command line and takes specific command line arguments

Input
#####
    * -i : A path to the ICGC .vcf file
    * -s : A schema file containing the field names in the annotations and to use for the output file
    * -o : A path to the output folder, where the transformed CSV data will go

Output
######
    * A .csv file with mutation data where each row contains one mutation and one unique annotation

Usage
#####
    * python convert_icgc_vcf.py -h

    *Gives a description of the neccessary commands

    * python convert_icgc_vcf.py -i <path/input_file.vcf> -s <path/schema.json> -o <path/>

    *Runs the script with the given input vcf and schema json and outputs a csv file


*Run map_icgc.py*
^^^^^^^^^^^^^^^^^

*Summary*
"""""""""

The python script map_icgc.py will take the output of the vcf convertor script and:
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

The script must be called from the command line and takes specific command line arguments

Input
#####
    * -i : A path to the ICGC .csv file
    * -m : A path to the folder containing mapping files
    * -d : The name of the doid mapping file
    * -e : The name of the ensp to uniprot accession mapping file
    * -o : A path to the output folder

Output
######
    * A .csv file with mutation data formatted to the biomuta field structure

Usage
#####
    * python map_icgc.py -h

    *Gives a description of the neccessary commands

    * python map_icgc.py -i <path/input_file.vcf> -m <path/> -d doid_mapping_file.csv -e enst_mapping_file.csv -o <path/>

    *Runs the script with the given csv file and outputs a csv file formatted for the final biomuta master file

*Additional Notes*
------------------

**Mapping**
###########

All the mapping files are alable in the scripts repository in the folder pipeline/convert_step2/mapping

The mapping files used for converting the ICGC csv are:

**DOID:** tcga_doid_mapping.csv

ICGC uses TCGA study terms, so the same TCGA to DOID parent terms are used for mapping (generated from previous Biomuta mapping):

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

Transcript ID (starts with ENST) was mapped to uniprot annotation accession

Mapping was NOT performed to uniprot canonical accession as this resulted in an issue with the final dataset in which a mutation for the same canonical accession would be listed with different amino acid changes