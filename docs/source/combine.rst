**Step 3: Combine**
===================

All resources are combined together into a master dataset. 

**Scripts**
-----------

combine_csv.py

**Procedure**
-------------

*Run combine_csv.py*
^^^^^^^^^^^^^^^^^^^^

*Summary*
"""""""""

All of the mutation data for each source was converted to a standardized data structure in the convert step.

Now, all of these separate csv files (one for each source) will be combined into a master csv

All csv files to be combined should be in a folder together with no additional csv files

*Script Specifications*
"""""""""""""""""""""""

The script must be called from the command line and takes specific command line arguments

Input
#####
    * -i : The folder containing csv mutation files to combine
    * -o : The folder to output the combined mutation file

Output
######
    * A csv file combining all csv files in a given folder

Usage
#####
    * python combine_csv.py -h

    *Gives a description of the neccessary commands

    * python combine_csv.py -i <path/> -o <path/>

    *Runs the script with the given folder and combines all csv files in that folder

*Additional Notes*
------------------

**Final fields**
^^^^^^^^^^^^^^^^

Field,Description
sample_name,Sample ID provided by the original resource (for v-5.0 only applies to TCGA and COSMIC)
chr_id,Chromosome number only (no 'chr')
start_pos,Genomic coordinates (For v-5.0 these are all with ref GRCh38)
end_pos,Identical to the start positoon because all mutations are Specifications
ref_nt,Reference nucleotide
alt_nt,Nucleotide mutation
aa_pos,Amino acide number of the amino acide change in the human_protein_transcriptlocus
ref_aa,Reference amino acid
alt_aa,Amino acid variation caused by the mutation
do_name,DO parent term
uniprot_canonical_ac,Uniprot accession for the specific ENST or ENSP listed from the source
source,Original data source of the mutation

