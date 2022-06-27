import os,sys
import string
import commands
import glob
import json
import csv
from optparse import OptionParser
import bz2


__status__ = "Dev"
__version__ = "4.0"






##############################################
def main():

    in_file_one = "/data/projects/biomuta/downloads/v-1.0/dbSNP/SNPdb_s3.tsv"
    in_file_two = "/data/projects/biomuta/downloads/v-1.0/dbSNP/2019_germline_mutations.csv"
            

    seen= {"snpaa":{}, "snpnt":{}}
    snp_nt_total, snp_aa_total = 0, 0
    with open(in_file_one, 'r') as f_in:
        for line in f_in:
            row = line.split("\t")
            snp_nt = "%s:%s:%s:%s" % (row[1], row[2], row[6], row[7])
            snp_aa = "%s:%s:%s:%s" % (row[4], row[8], row[9], row[10])
            if snp_nt not in seen["snpnt"]:
                seen["snpnt"][snp_nt] = True
                snp_nt_total += 1
            if snp_aa not in seen["snpaa"]:
                seen["snpnt"][snp_aa] = True
                snp_aa_total += 1
    print snp_nt_total, snp_aa_total

    snp_aa_total = 0
    with open(in_file_two, 'r') as FR:
        row_count = 0
        for line in FR:
            row_count += 1
            row = line.strip().split(",")
            snp_aa = "%s:%s:%s:%s" % (row[1], row[2], row[3], row[4])
            if snp_aa not in seen["snpaa"]:
                seen["snpnt"][snp_aa] = True
                snp_aa_total += 1
    print snp_aa_total




if __name__ == '__main__':
	main()


