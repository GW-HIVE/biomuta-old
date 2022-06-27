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


    seen= {"uniprot2refseq":{}, "refseq2uniprot":{}, "refseq2ntsnp":{}, "rowmap":{}}
    in_file_one = "/data/projects/biomuta/downloads/v-1.0/dbSNP/SNPdb_s3.tsv"
    in_file_two = "/data/projects/biomuta/downloads/v-1.0/dbSNP/2019_germline_mutations.csv"
    with open(in_file_one, 'r') as FR:
        row_count = 0
        for line in FR:
            row_count += 1
            row = line.strip().split("\t")
            if row[0][0:1] == "#":
                continue
            snp_nt = "%s:%s:%s:%s" % (row[1], row[2], row[6], row[7])
            
            snp_aa_uniprot = "%s:%s" % (row[-2], row[-1])
            snp_aa_refseq = "%s:%s" % (row[4], row[8])
            #snp_aa_uniprot = "%s" % (row[-2])
            #snp_aa_refseq = "%s" % (row[4])
            
            seen["uniprot2refseq"][snp_aa_uniprot] = snp_aa_refseq
            seen["refseq2uniprot"][snp_aa_refseq] = snp_aa_uniprot
            seen["refseq2ntsnp"][snp_aa_refseq] = snp_nt

    row_count = 0
    n1, n2 = 0, 0
    with open(in_file_two, 'r') as FR:
        csv_grid = csv.reader(FR, delimiter=',', quotechar='"')
        for row in csv_grid:
            row_count += 1
            if row_count == 1:
                newrow = row + ["refseq_ac", "refseq_pos",  "pos_shift","map_flag", "rsid_list"]
                print "\"%s\"" % ( "\",\"".join(newrow))
                continue
            snp_aa_uniprot = "%s:%s" % (row[1], row[2])
            #snp_aa_uniprot = "%s" % (row[1])
            
            if snp_aa_uniprot not in seen["rowmap"]:
                seen["rowmap"][snp_aa_uniprot] = []
            seen["rowmap"][snp_aa_uniprot].append(row)

            if snp_aa_uniprot in seen["uniprot2refseq"]:
                snp_aa_refseq = seen["uniprot2refseq"][snp_aa_uniprot]
                n1 += 1
            else:
                n2 += 1

    snpaa2rsid = {}
    for in_file in glob.glob("outdir/refsnp-chr*.csv"):
        with open(in_file, 'r') as f_in:
            for line in f_in:
                row = line.split(",")
                snp_nt = ":".join(row[2:6])
                snp_aa_refseq = ":".join(row[-4:-2])
                #snp_aa_refseq = ":".join(row[-4:-3])
                if snp_aa_refseq not in snpaa2rsid:
                    snpaa2rsid[snp_aa_refseq] = []
                snpaa2rsid[snp_aa_refseq].append(row[0])


    n1, n2 = 0, 0 
    for snp_aa_uniprot in seen["rowmap"]:
        snp_aa_refseq_one = seen["uniprot2refseq"][snp_aa_uniprot] if snp_aa_uniprot in seen["uniprot2refseq"] else ""
        refseq_ac,refseq_pos_one = snp_aa_refseq_one.split(":")[0], snp_aa_refseq_one.split(":")[1]
        #refseq_ac,refseq_pos_one = snp_aa_refseq_one, 0
        
        refseq_pos_two = int(refseq_pos_one) - 1
        snp_aa_refseq_two = refseq_ac + ":" + str(refseq_pos_two)
        
        refseq_pos_three = int(refseq_pos_one) + 1
        snp_aa_refseq_three = refseq_ac + ":" + str(refseq_pos_three)

        rsid_list = []
        refseq_pos = refseq_pos_one
        flag = False
        shift = ""
        if snp_aa_refseq_one in snpaa2rsid:
            rsid_list = snpaa2rsid[snp_aa_refseq_one]
            flag = True
            shift = "0"
        elif snp_aa_refseq_two in snpaa2rsid:
            rsid_list = snpaa2rsid[snp_aa_refseq_two]
            refseq_pos = refseq_pos_two
            flag = True
            shift = "-1"
        elif snp_aa_refseq_three in snpaa2rsid:
            rsid_list = snpaa2rsid[snp_aa_refseq_three]
            refseq_pos = refseq_pos_three
            flag = True
            shift = "+1"

        for row in seen["rowmap"][snp_aa_uniprot]:
            newrow = row + [refseq_ac,str(refseq_pos),shift,str(flag), "|".join(rsid_list)]
            print "\"%s\"" % ( "\",\"".join(newrow))



if __name__ == '__main__':
	main()


