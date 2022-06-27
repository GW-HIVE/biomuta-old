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


    seen= {"snpaa":{}, "snpnt":{}}
    snp_nt_total, snp_aa_total = 0, 0
    for in_file in glob.glob("outdir/refsnp-chr*.csv"):
        snp_nt_count, snp_aa_count = 0, 0
        with open(in_file, 'r') as f_in:
            for line in f_in:
                row = line.split(",")
                snp_nt = ":".join(row[2:6])
                snp_aa = ":".join(row[-4:])
                if snp_nt not in seen["snpnt"]:
                    seen["snpnt"][snp_nt] = True
                    snp_nt_count += 1
                if snp_aa not in seen["snpaa"]:
                    seen["snpnt"][snp_aa] = True
                    snp_aa_count += 1
        snp_nt_total += snp_nt_count
        snp_aa_total += snp_aa_count
        print in_file, snp_nt_count, snp_aa_count
    print "Total: ", snp_nt_count, snp_aa_count




if __name__ == '__main__':
	main()


