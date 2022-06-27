import os,sys
import string
import csv
import json
from optparse import OptionParser
from Bio import SeqIO
from Bio.Seq import Seq


__version__="1.0"
__status__ = "Dev"


###############################
def translate_codon_one(codon_seq, codon_dic, last_codon_flag):
    aa = ""
    if codon_seq in codon_dic:
        if last_codon_flag and codon_dic[codon_seq] == "*":
            aa = ""
        elif codon_dic[codon_seq] == "*":
            aa = "*"
        else:
            aa = codon_dic[codon_seq]
    elif "N" in codon_seq:
        aa = "X"
	
    return aa

###############################
def translate_codon_two(codon_seq, codon_dic):
    aa = ""
    if codon_seq in codon_dic:
        aa = codon_dic[codon_seq]
    elif "N" in codon_seq:
        aa = "X"

    return aa




####################################
def load_mut_list_in_cds(mutlist_file, chr_list, cdslocus2tridlist):

        mut_list = []
	mut_count_hash = {"total":0, "coding":0, "noncoding":0, "incontigs":0}
        with open(mutlist_file, 'rb') as tsvfile:
                tsvreader = csv.reader(tsvfile, delimiter=',', quotechar='"')
                row_count = 0
                for row in tsvreader:
                        row_count += 1
                        if row_count == 1:
                                continue
			mut_count_hash["total"] += 1
                        if row[0] not in chr_list:
				mut_count_hash["incontigs"] += 1
                                continue
                        if len(row[3]) > 1 or len(row[4]) > 1:
				continue
			k = row[0] + ":" + row[1]
                        if k not in cdslocus2tridlist:
                                mut_count_hash["noncoding"] += 1
				continue
			val_list = [str(row_count), row[0], row[1],row[3], row[4]]
                        mut_list.append(val_list)
			mut_count_hash["coding"] += 1
        return mut_list, mut_count_hash

####################################
def get_codon_info(exon_objs, trs_id, exon_index, mut_pos_in_exon, mut_pos_in_codon):


	obj = exon_objs[trs_id]["exonlist"][exon_index]
	exon_len = obj["stop"] - obj["start"] + 1
	prev_obj, next_obj = {"seq":"NNNN"}, {"seq":"NNNN"}
	if exon_index > 0:
		prev_obj = exon_objs[trs_id]["exonlist"][exon_index-1]
	if exon_index < len(exon_objs[trs_id]["exonlist"]) - 1:
		next_obj = exon_objs[trs_id]["exonlist"][exon_index+1]

	codon_start_in_exon = mut_pos_in_exon - mut_pos_in_codon
	codon_stop_in_exon = codon_start_in_exon + 3
	codon_start_in_cds = obj["cdsstart"] - 1 + codon_start_in_exon
	codon_stop_in_cds = codon_start_in_cds + 3

	codon_seq = ""
	if codon_start_in_exon < 0:
		codon_seq = prev_obj["seq"][codon_start_in_exon:]
		codon_seq += obj["seq"][0:obj["frame"]]
		codon_start_in_cds += codon_start_in_exon
		codon_stop_in_cds = codon_start_in_cds + 3
	elif codon_start_in_exon  ==  exon_len - 2:
		codon_seq = obj["seq"][codon_start_in_exon:] + next_obj["seq"][0:1]
	elif codon_start_in_exon  ==  exon_len - 1:
		codon_seq = obj["seq"][codon_start_in_exon:] + next_obj["seq"][0:2]
	else:
		codon_seq = obj["seq"][codon_start_in_exon:codon_stop_in_exon]

	return codon_seq, codon_start_in_cds






###############################
def main():

            
    dna_seq_file = "/data/external/ucsc/grch37/targets/Homo_sapiens.GRCh37.75.dna.all.fa"
    cds_seq_file = "/data/external/ucsc/grch37/targets/Homo_sapiens.GRCh37.75.cds.all.fa"
    pep_seq_file = "/data/external/ucsc/grch37/targets/Homo_sapiens.GRCh37.75.pep.all.fa"
    codon_file = "/data/external/misc/codon_list.txt"
    exon_list_file = "outdir/exonlist.chr12.json"

    chr_list = ["12"]

    #Load cds sequences
    cds_seq_hash = {}
    for record in SeqIO.parse(cds_seq_file, "fasta"):
        chr_id = record.description.split(" ")[2].split(":")[2]
        frame = record.description.split(" ")[2].split(":")[-1]
        if chr_id not in chr_list:
            continue
        trs_id = record.id.split(".")[0]
        cds_seq_hash[trs_id] = record.seq.upper()

    #Load pep sequences
    pep_seq_hash = {}
    trsid2pepid = {}
    for record in SeqIO.parse(pep_seq_file, "fasta"):
        desc = record.description
        trs_id = desc.split("transcript:")[1].split()[0].split(".")[0]
        pep_seq_hash[trs_id] = str(record.seq.upper())
        trsid2pepid[trs_id] =  record.id.split(".")[0]


    genome_seq_hash = {}
    for record in SeqIO.parse(dna_seq_file, "fasta"):
        chr_id = record.id.strip().replace("chr", "")
        if chr_id in chr_list:
            genome_seq_hash[chr_id] = record.seq.upper()


    #Load codon dict
    codon_dic = {}
    with open(codon_file, "r") as tsvfile:
        for row in tsvfile:
            row = row.strip().split('\t')
            codon_dic[row[0]] = row[1]

    exon_objs = json.loads(open(exon_list_file, "r").read())

    #Populate pos to trsid mapping
    cdslocus2tridlist = {}
    for trs_id in exon_objs:
        if exon_objs[trs_id]["exonlist"][0]["chr_id"] not in chr_list:
            continue
        for obj in exon_objs[trs_id]["exonlist"]:
            for j in xrange(int(obj["start"]), int((obj["stop"])) + 1):
                k = obj["chr_id"] + ":" + str(j)
                if k not in cdslocus2tridlist:
                    cdslocus2tridlist[k] = []
                if trs_id not in cdslocus2tridlist[k]:
                    cdslocus2tridlist[k].append(trs_id)


    mut = "12:25398284:C>A"

    chr_id = mut.split(":")[0]
    chr_pos = int(mut.split(":")[1])
    ref = mut.split(":")[2].split(">")[0]
    alt = mut.split(":")[2].split(">")[1]

    cmp_hash = {"A":"T", "T":"A", "C":"G", "G":"C"}



    print "in_genome,in_cds,in_peptide"
    locus = chr_id + ":" + str(chr_pos)

    for trs_id in cdslocus2tridlist[locus]:
        frame1 = exon_objs[trs_id]["exonlist"][0]["frame"]
        frame2 = exon_objs[trs_id]["exonlist"][-1]["frame"]
        for exon_index in xrange(0, len(exon_objs[trs_id]["exonlist"])):
            obj = exon_objs[trs_id]["exonlist"][exon_index]
            alt_new = cmp_hash[alt] if obj["strand"] == "-" else alt
            ref_new = cmp_hash[ref] if obj["strand"] == "-" else ref
            if obj["status"] == "removed":
                continue
            if chr_pos < int(obj["start"]) or chr_pos > int(obj["stop"]):
                continue
            exon_len = obj["stop"] - obj["start"] + 1
            mut_pos_in_exon = chr_pos - obj["start"]
            if obj["strand"] == "-":
                mut_pos_in_exon = obj["stop"] - chr_pos
            mut_pos_in_cds = obj["cdsstart"] - 1 + mut_pos_in_exon
            mut_pos_in_codon = (mut_pos_in_exon-obj["frame"])%3
            ref_in_genome = genome_seq_hash[chr_id][chr_pos-1]
            ref_in_exon = obj["seq"][mut_pos_in_exon]
            ref_in_cds = cds_seq_hash[trs_id][mut_pos_in_cds]
            codon_index_in_cds = int(mut_pos_in_cds/3)
            ind = codon_index_in_cds*3
            codon_seq = str(cds_seq_hash[trs_id][ind:ind+3])
            last_codon_flag = codon_index_in_cds == len(pep_seq_hash[trs_id]) 
            aa_translated = translate_codon_two(codon_seq, codon_dic)
            if aa_translated == "":
                continue
            aa_in_pep_seq = pep_seq_hash[trs_id][codon_index_in_cds] if last_codon_flag == False else "*"
            if aa_translated != aa_in_pep_seq:
                continue
            alt_codon_seq = codon_seq[:mut_pos_in_codon]
            alt_codon_seq += alt_new + codon_seq[mut_pos_in_codon+1:]
            aa_alt = translate_codon_two(alt_codon_seq, codon_dic)
            out_row = ["%s:%s:%s>%s" % (chr_id,chr_pos,ref,alt)]
            out_row += ["%s:%s:%s>%s" % (trs_id,mut_pos_in_cds + 1,codon_seq,alt_codon_seq)]
            out_row += ["%s:%s:%s>%s" % (trsid2pepid[trs_id],codon_index_in_cds + 1,aa_translated,aa_alt)]
            print ",".join(out_row)

                







if __name__ == '__main__':
    main()








