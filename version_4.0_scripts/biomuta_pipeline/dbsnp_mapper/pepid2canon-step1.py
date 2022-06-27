import os,sys
import string
import csv
import json
from Bio import SeqIO
from Bio.Seq import Seq
import commands
import subprocess
from optparse import OptionParser


__version__="1.0"
__status__ = "Dev"





###############################
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
    uniprot_ver = config_json["uniprotversion"]


    
    pep_seq_file = config_json["ref"][grch_ver]["pepseqfile"]
    canon_seq_file = config_json["ref"][uniprot_ver]["canonseqfile"]
    genelocus_file = config_json["ref"][uniprot_ver]["genelocusfile"]


    seq_hash = {}
    for record in SeqIO.parse(pep_seq_file, "fasta"):
        seq_id = record.id.split(".")[0]
        seq_hash[seq_id] = str(record.seq.upper())


    for record in SeqIO.parse(canon_seq_file, "fasta"):
        seq_id = record.id.split("|")[1]
        seq_hash[seq_id] = str(record.seq.upper())

    FW = open("outdir/pepid2isoformac.csv", "w")
    FW.write("pep_id,flag,uniprotkb_canonical_ac\n")
    with open(genelocus_file, "r") as FR:
        data_grid = csv.reader(FR, delimiter=',', quotechar='"')
        row_count = 0
        for row in data_grid:
            row_count += 1
            if row_count > 1:
                canon, isoform, pep_id = row[0], row[1], row[3]
                if canon not in seq_hash:
                    print "No sequence found for %s" % (canon)
                    continue
                if pep_id not in seq_hash:
                    print "No sequence found for %s" % (seq_id)
                    continue

                if canon == isoform:
                    flag = "exact" if seq_hash[canon] == seq_hash[pep_id] else "nonexact"
                    FW.write("%s,%s,%s\n" % (pep_id, flag, canon))
    FW.close()





        
        
                        


if __name__ == '__main__':
        main()



