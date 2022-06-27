import os,sys
import string
import commands
import glob
import json
import csv
from optparse import OptionParser


__status__ = "Dev"
__version__ = "4.0"






##############################################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-i","--configfile",action="store",dest="configfile",help="NT file")
    
    (options,args) = parser.parse_args()
    for file in ([options.configfile]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    config_json = json.loads(open(options.configfile, "r").read())
    grch_ver = config_json["grchversion"]
    vcf_file = config_json["ref"][grch_ver]["dbsnpfile"]


    chr_list = map(str,range(1,23)) + ["X","Y"]
    #chr_list = ["1"]
    chrac2id = {"NC_000023":"X", "NC_000024":"Y"}
    for chr_id in map(str,range(1,23)) :
        ac = "NC_00000" + chr_id if len(chr_id) == 1 else "NC_0000" + chr_id
        chrac2id[ac] = chr_id 

    key_list = ["CAF","SAO","GENEINFO","VC", "NSF", "NSM", "NSN", "SYN"]

    FW = open("outdir/dbsnp.1.csv", "w")
    field_list = ["chr", "pos", "id", "ref", "alt", "caf", "sao", "geneinfo","vc", "nsf","nsm", "nsn", "syn"]
    FW.write("\"%s\"\n" % ("\",\"".join(field_list)))
    with open(vcf_file, 'r') as FR:
        row_count = 0
        for line in FR:
            row_count += 1
            row = line.strip().split("\t")
            if row[0][0:1] == "#":
                continue
            chr_id = row[0]
            chr_ac = row[0].split(".")[0]
            chr_id = chrac2id[chr_ac] if chr_ac in chrac2id else row[0]

            tv_list = []
            tv_list.append(chr_id in chr_list)
            tv_list.append(len(row[3]) == 1 and len(row[4]) == 1)
            tv_list.append(row[3] in ["A","C", "G", "T"] and row[4] in ["A","C", "G", "T"])
            if False in tv_list:
                continue

	    val_list = [chr_id,row[1],row[2],row[3],row[4]]
            mut_id = row[2]
            if mut_id not in ["", "."]:
                val_dict = {}
                for s in row[7].split(";"):
                    k, v = "", ""
		    if s.find("=") != -1:
                        k, v = s.split("=")[0], "=".join(s.split("=")[1:])
                    else:
			k, v = s.strip(), s.strip()
		    val_dict[k] = v
                for k in key_list:
                   val_list.append(val_dict[k] if k in val_dict else "")
                FW.write("\"%s\"\n" % ("\",\"".join(val_list)))
    FW.close()







if __name__ == '__main__':
	main()


