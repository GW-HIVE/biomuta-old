import os,sys
import string
import csv
import json
from Bio import SeqIO
from Bio.Seq import Seq
import commands
import subprocess


__version__="1.0"
__status__ = "Dev"





###############################
def main():

    remote_dir = "http://dev.data.glygen.org/datasets/unreviewed/"
    local_dir = "/data/external/uniprot/2018_11/"
    
    #remote_dir = "http://dev.data.glygen.org/datasets/reviewed/"
    #local_dir = "/data/external/uniprot/2017_11/"


    file_list = ["human_protein_canonicalsequences.fasta", "human_protein_allsequences.fasta"]
    for file_name in file_list:
        remote_path = "%s/%s" % (remote_dir, file_name)
        local_path = "%s/targets/%s" % (local_dir, file_name)
        cmd = "/usr/bin/wget %s -O %s " % (remote_path, local_path) 
        x = commands.getoutput(cmd)

    file_list = ["human_protein_genelocus.csv", "human_protein_idmapping.csv"]
    for file_name in file_list:
        remote_path = "http://dev.data.glygen.org/datasets/unreviewed/%s" % (file_name)
        local_path = "%s/idmap/%s" % (local_dir, file_name)
        cmd = "/usr/bin/wget %s -O %s " % (remote_path, local_path)
        x = commands.getoutput(cmd)








        
        
                        


if __name__ == '__main__':
        main()



