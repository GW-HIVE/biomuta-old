# ===========================================================================
#
#                            PUBLIC DOMAIN NOTICE
#               National Center for Biotechnology Information
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the author's official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software/database is freely available
#  to the public for use. The National Library of Medicine and the U.S.
#  Government have not placed any restriction on its use or reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, the NLM and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. The NLM and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any particular
#  purpose.
#
#  Please cite the author in any work or product based on this material.
#
# ===========================================================================
# Script name: rsjson_demo.py
# Description: a demo script to parse dbSNP RS JSON object.  The script will
# produce tab-delimited output containing tthe assembly version, sequence ID,
# position, reference allele, variant allele and ClinVar clinical significance
# if available.
# Author:  Lon Phan  lonphan@ncbi.nlm.nih.gov
# For help please contact: tkt-varhd@ncbi.nlm.nih.gov
#
#
# ---------------------------------------------------------------------------


import argparse
import json
import bz2


def printAllele_annotations(primary_refsnp):
    '''
    rs clinical significance
    '''
    for annot in primary_refsnp['allele_annotations']:
        for clininfo in annot['clinical']:
            print(",".join(clininfo['clinical_significances']))


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



parser = argparse.ArgumentParser(description='Example of parsing JSON RefSNP Data')
parser.add_argument('-i', dest='input_fn', required=True, help='The name of the input file to parse')

args = parser.parse_args()


cnt = 0
with bz2.BZ2File(args.input_fn, 'rb') as f_in:
    for line in f_in:
        rs_obj = json.loads(line.decode('utf-8'))
        #print(rs_obj['refsnp_id'] + "\t")  # rs ID
        rs_id = rs_obj['refsnp_id']
        if 'primary_snapshot_data' in rs_obj:
            rows = printPlacements(rs_obj['primary_snapshot_data']['placements_with_allele'])
            #printAllele_annotations(rs_obj['primary_snapshot_data'])
            for row in rows:
                for o_one in rs_obj["primary_snapshot_data"]["allele_annotations"]:
                    for o_two in o_one["assembly_annotation"]:
                        for o_three in o_two["genes"]:
                            for o_four in o_three["rnas"]:
                                if "protein" in o_four:
                                    v = o_four["protein"]["variant"]["spdi"]
                                    refseq_ac, pos = v["seq_id"], v["position"]
                                    ref_aa, alt_aa = v["deleted_sequence"], v["inserted_sequence"]
                                    newrow = [rs_id] + row + [refseq_ac, pos, ref_aa, alt_aa]
                                    if ref_aa != alt_aa:
                                        print newrow

        cnt = cnt + 1
