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

        
    remote_dir = "ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/"
    local_dir = "/data/external/ucsc/grch38/genome/"

    chr_list = map(str,range(1,23)) + ["X","Y","MT"]
    chr_list = []
    for chr_id in chr_list:
        remote_path = "%s/Homo_sapiens.GRCh38.dna.chromosome.%s.fa.gz" % (remote_dir, chr_id)
        local_path = "%s/chromosome/Homo_sapiens.GRCh38.dna.chromosome.%s.fa.gz" % (local_dir, chr_id)
        cmd = "/usr/bin/wget %s -O %s " % (remote_path, local_path) 
        x = commands.getoutput(cmd)
        cmd = "/usr/bin/gunzip %s" % (local_path)
        x = commands.getoutput(cmd)

    
    cmd = "cat /data/external/ucsc/grch38/genome/chromosome/Homo_sapiens.GRCh38.dna.chromosome.*.fa > "
    cmd += "/data/external/ucsc/grch38/targets/Homo_sapiens.GRCh38.95.dna.all.fa "
    x = commands.getoutput(cmd)


        
        
                        


if __name__ == '__main__':
        main()



