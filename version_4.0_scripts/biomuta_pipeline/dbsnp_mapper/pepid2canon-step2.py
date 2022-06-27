import os,sys
import string
import csv
import json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
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

  
    seq_hash = {}
    for record in SeqIO.parse(pep_seq_file, "fasta"):
        seq_id = record.id.split(".")[0]
        seq_hash[seq_id] = str(record.seq.upper())

    for record in SeqIO.parse(canon_seq_file, "fasta"):
        seq_id = record.id.split("|")[1]
        seq_hash[seq_id] = str(record.seq.upper())

    FW = open("outdir/map2uniprot.ranges.csv", "w")
    FW.write("peptide_id,uniprotkb_canonical_ac,ranges\n")
    with open("outdir/pepid2isoformac.csv", "r") as FR:
        data_grid = csv.reader(FR, delimiter=',', quotechar='"')
        row_count = 0
        for row in data_grid:
            row_count += 1
            if row_count == 1:
                continue
            pep_id, flag, canon = row[0], row[1], row[2]
            if flag in ["exact", "ignore"]:
                continue
            seq1, seq2 = seq_hash[pep_id], seq_hash[canon]
            aln_list = pairwise2.align.globalxx(seq1, seq2)
            pos1, pos2 = 0, 0
            prev_pos1, prev_pos2 = 0, 0
            prev_aa1, prev_aa2 = "-", "-"
            start_pos1, start_pos2, end_pos1, end_pos2 = 0, 0, 0, 0
            range_list = []
            alnlen = 0
            for i in xrange(0, len(aln_list[0][0])):
                aa1, aa2 = aln_list[0][0][i], aln_list[0][1][i]
                pos1 += 1 if aa1 != "-" else 0
                pos2 += 1 if aa2 != "-" else 0
                c1 = aa1 != "-" and aa2 != "-"
                c2 = prev_aa1 == "-" or prev_aa2 == "-"
                if c1 and c2:
                    start_pos1, start_pos2 = pos1, pos2
                c1 = aa1 == "-" or aa2 == "-"
                c2 = prev_aa1 != "-" and prev_aa2 != "-"
                if (c1 and c2) or pos1 == len(seq1):
                    end_pos1, end_pos2 = prev_pos1, prev_pos2
                    alnlen += end_pos1 - start_pos1 + 1
                    r = "%s-%s:%s-%s"%(start_pos1,end_pos1,start_pos2,end_pos2)
                    range_list.append(r)
                prev_aa1 = aa1
                prev_aa2 = aa2
                prev_pos1 = pos1
                prev_pos2 = pos2
            ratio1 = float(alnlen)/len(seq1) 
            ratio2 = float(alnlen)/len(seq2)
            if ratio1 > 0.7 and len(range_list) <= 10:
                FW.write("%s,%s,%s\n" % (pep_id, canon, "|".join(range_list)))
    FW.close()







if __name__ == '__main__':
        main()





