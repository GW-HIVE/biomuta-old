**Convert: CIVIC**
==================

**Scripts**
-----------

**Procedure**
-------------

*Perform liftover of mutations from GRCh37 to GRCh38*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*Summary*
"""""""""

The most recent data release for CIVIC is aligned to the GRCH37 human reference genome. For this update however, we are using the human reference genome GRCh38.

To convert coordinates between the two reference genomes, we use a 'liftover' tool to remap the genomic coordinates. 

The CIVIC file is very small in size, so we can use the ENSEMBL online liftover tool: https://useast.ensembl.org/Homo_sapiens/Tools/AssemblyConverter?db=core

Run the downloaded VCF throught the tool with the default parameters (except chnaging the file type to VCF)

Redownload the transformed VCF and use that VCF for the next step. 

*Run convert_civic_vcf.py*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*Summary*
"""""""""

The python scipt convert_civic_vcf.py will convert the vcf formatted file to a csv file. 

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

**Running convert_civic_vcf.py**
################################

The script must be called from the command line and takes specific command line arguments



