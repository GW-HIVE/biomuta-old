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


    local_dir = "/data/external/ucsc/hg38/targets/"

    for seq_type in ["cds", "pep"]:
        remote_dir = "ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/%s/" % (seq_type)
        remote_path = "%s/Homo_sapiens.GRCh38.%s.all.fa.gz" % (remote_dir, seq_type)
        local_path = "%s/Homo_sapiens.GRCh38.%s.all.fa.gz" % (local_dir, seq_type)
        cmd = "/usr/bin/wget %s -O %s " % (remote_path, local_path) 
        x = commands.getoutput(cmd)
        cmd = "/usr/bin/gunzip %s" % (local_path)
        x = commands.getoutput(cmd)


    remote_path = "ftp://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.95.gtf.gz"
    local_path = "/data/external/ucsc/hg38/gtf/Homo_sapiens.GRCh38.95.gtf.gz"
    cmd = "/usr/bin/wget %s -O %s " % (remote_path, local_path)
    x = commands.getoutput(cmd)
    cmd = "/usr/bin/gunzip %s" % (local_path)
    x = commands.getoutput(cmd)







        
        
                        


if __name__ == '__main__':
        main()



