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



def printPlacements(info):
    '''
    rs genomic positions
    '''

    rows = []
    for alleleinfo in info:
        # has top level placement (ptlp) and assembly info
        if alleleinfo['is_ptlp'] and \
                len(alleleinfo['placement_annot']['seq_id_traits_by_assembly']) > 0:
            assembly_name = alleleinfo['placement_annot'] \
                                      ['seq_id_traits_by_assembly'] \
                                      [0]['assembly_name']

            for a in alleleinfo['alleles']:
                spdi = a['allele']['spdi']
                if spdi['inserted_sequence'] != spdi['deleted_sequence']:
                    (ref, alt, pos, seq_id) = (spdi['deleted_sequence'],
                                               spdi['inserted_sequence'],
                                               spdi['position'],
                                               spdi['seq_id'])
                    rows.append([assembly_name, seq_id, str(pos), ref, alt])
                    break
    return rows






##############################################
def main():


    seen= {"aasnpmap":{}, "rowmap":{}, "snpaa2snpnt":{}}

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
            snp_aa_uniprot = "%s:%s:%s:%s" % (row[-2], row[-1], row[9], row[10])
            snp_aa_refseq = "%s:%s:%s:%s" % (row[4], row[8], row[9], row[10])
            seen["aasnpmap"][snp_aa_uniprot] = snp_aa_refseq
            seen["aasnpmap"][snp_aa_refseq] = snp_aa_uniprot
            seen["snpaa2snpnt"][snp_aa_refseq] = snp_nt


    with open(in_file_two, 'r') as FR:
        row_count = 0
        for line in FR:
            row_count += 1
            row = line.strip().split(",")
            if row_count == 1:
                print ",".join(row + ["SNP", "RSID"])
                continue
            snp_aa_uniprot = "%s:%s:%s:%s" % (row[1], row[2], row[3], row[4])
            if snp_aa_uniprot in seen["aasnpmap"]:
                snp_aa_refseq = seen["aasnpmap"][snp_aa_uniprot]
                seen["rowmap"][snp_aa_refseq] = row

    #for snp_aa in seen["aasnptwo"]:
    #    print snp_aa, snp_aa in seen["aasnpmap"]
    #sys.exit()

    FW = open("outdir/JUNK.csv", "w")
    for in_file in glob.glob("/data/projects/biomuta/downloads/v-1.0/dbSNP/refsnp-chr*.json.bz2"):
        with bz2.BZ2File(in_file, 'rb') as f_in:
            for line in f_in:
                rs_obj = json.loads(line.decode('utf-8'))
                rs_id = rs_obj['refsnp_id']
                if 'primary_snapshot_data' in rs_obj:
                    rows = printPlacements(rs_obj['primary_snapshot_data']['placements_with_allele'])
                    for row in rows:
                        for o_one in rs_obj["primary_snapshot_data"]["allele_annotations"]:
                            for o_two in o_one["assembly_annotation"]:
                                for o_three in o_two["genes"]:
                                    for o_four in o_three["rnas"]:
                                        if "protein" in o_four:
                                            if "spdi" not in o_four["protein"]["variant"]:
                                                continue
                                            v = o_four["protein"]["variant"]["spdi"]
                                            refseq_ac, pos = v["seq_id"].split(".")[0], v["position"]
                                            ref_aa, alt_aa = v["deleted_sequence"], v["inserted_sequence"]
                                            newrow = [str(rs_id)] + row + [refseq_ac, str(pos), ref_aa, alt_aa]
                                            if ref_aa != alt_aa:
                                                snp_aa_refseq = "%s:%s:%s:%s" % (refseq_ac,pos,ref_aa,alt_aa)
                                                #print snp_aa_refseq, snp_aa_refseq in seen["aasnpmap"], newrow
                                                if snp_aa_refseq in seen["aasnpmap"]:
                                                    FW.write("%s\n" % (",".join(newrow)))
                                                    #print snp_aa_refseq, snp_aa_refseq in seen["rowmap"]
                                                    #if snp_aa_refseq in seen["rowmap"]:
                                                    #    snp_nt = seen["snpaa2snpnt"][snp_aa_refseq]
                                                    #    print seen["rowmap"][snp_aa_refseq] + [rs_id, snp_nt]

    FW.close()

if __name__ == '__main__':
	main()


