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

    FW = open("outdir/dbsnp-aansp.csv", "w")
    for in_file in glob.glob("/data/projects/biomuta/downloads/v-1.0/dbSNP/refsnp-chr*.json.bz2"):
        file_name = in_file.split("/")[-1].split(".")[0]
        out_file = "outdir/%s.csv" %(file_name)
        print "Now parsing ", in_file 
        FW = open(out_file, "w")
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
                                            FW.write("%s\n" % (",".join(newrow)))
        FW.close()

if __name__ == '__main__':
	main()


