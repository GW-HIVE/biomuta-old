import os,sys
import string
import csv
import json
from optparse import OptionParser
from Bio import SeqIO
from Bio.Seq import Seq


__version__="1.0"
__status__ = "Dev"



def dump_formatted_cds (trs_id, pep_id, cds_seq, pep_seq):

    cds_pos = 0
    pep_pos = 0
    buffer = ">%s (%s)" % (trs_id, pep_id)
    while cds_pos < len(cds_seq) - 2:
        codon = cds_seq[cds_pos:cds_pos+3]
        aa = pep_seq[pep_pos] if pep_pos < len(pep_seq) else "*"
        nl = "\n" if pep_pos  % 10 == 0 else ""
        buffer += "%s%s(%s) " % (nl,codon,aa)
        cds_pos += 3
        pep_pos += 1
   
    return buffer.strip()




###############################
def main():


    cds_seq_file = "/data/external/ucsc/grch37/targets/Homo_sapiens.GRCh37.75.cds.all.fa"
    pep_seq_file = "/data/external/ucsc/grch37/targets/Homo_sapiens.GRCh37.75.pep.all.fa"

    #Load cds sequences
    cds_seq_hash = {}
    for record in SeqIO.parse(cds_seq_file, "fasta"):
        trs_id = record.id.split(".")[0]
        cds_seq_hash[trs_id] = record.seq.upper()

    pep_seq_hash = {}
    trsid2pepid = {}
    for record in SeqIO.parse(pep_seq_file, "fasta"):
        desc = record.description
        trs_id = desc.split("transcript:")[1].split()[0].split(".")[0]
        pep_seq_hash[trs_id] = str(record.seq.upper())
        trsid2pepid[trs_id] =  record.id.split(".")[0]

    trs_id = "ENST00000311936"

    print dump_formatted_cds (trs_id, trsid2pepid[trs_id], cds_seq_hash[trs_id], pep_seq_hash[trs_id])
                







if __name__ == '__main__':
    main()








