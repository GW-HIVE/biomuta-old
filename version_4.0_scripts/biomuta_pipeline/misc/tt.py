import os,sys
import string
import commands
from optparse import OptionParser
import glob
import json
import csv


__status__ = "Dev"
__version__ = "4.0"






##############################################
def main():

    chr_list = ["12"]

    src_type = "tcga"
    input_dir = "/data/projects/biomuta/downloads/v-4.0/tcga/"
    pattern = input_dir + "/*/*/*.vcf"
    #pattern = "/data/projects/biomuta/downloads/v-4.0/tcga/ACC/a9*/*.vcf" 
    file_list = glob.glob(pattern)

    start_pos = 25357723 
    end_pos = 25403870

    file_count = 1
    for vcf_file in file_list:
        cancer_type = vcf_file.split("/")[-3]
        patient_id = vcf_file.split("/")[-1].split(".")[0]
        out_file =  "tmp/%s-subject-%s.vcf" % (cancer_type, file_count)
        FW = open(out_file, "w")
        file_count += 1
        with open(vcf_file, 'r') as FR:
            row_count = 0
            for line in FR:
                row_count += 1
                row = line.strip().split("\t")
                chr_id = row[0]
                pos = int(row[1]) if line[0] != "#" else 0
                if line[0] == "#":
                    FW.write("%s\n" % (line.strip()))
                elif chr_id in chr_list and pos >= start_pos and pos <= end_pos:
                    FW.write("%s\n" % (line.strip()))

        FW.close()
	



if __name__ == '__main__':
	main()


