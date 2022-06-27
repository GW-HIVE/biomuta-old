import os,sys
import string
import csv
import json
from optparse import Option_parser
from Bio import Seq_iO
from Bio.Seq import Seq


__version__="1.0"
__status__ = "Dev"

###############################
def translate_codon(codon_seq, codon_dic, last_codon_flag):

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



#################################
def report_progress(progress_file, msg, open_new=False):

	msg = msg.strip()
        if open_new == True:
                with open(progress_file, "w") as FW:
                        FW.write(msg + "\n")
        else:
                with open(progress_file, "a") as FA:
                        FA.write(msg + "\n")
        return


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
                        if row[1] not in chr_list:
				mut_count_hash["incontigs"] += 1
                                continue
                        k = row[1] + ":" + row[2]
                        if k not in cdslocus2tridlist:
                                mut_count_hash["noncoding"] += 1
				continue
                        mut_list.append(row)
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

	usage = "\n%prog  [options]"
	parser = Option_parser(usage,version="%prog " + __version__)
	parser.add_option("-i","--configfile",action="store",dest="configfile",help="NT file")


	(options,args) = parser.parse_args()
	for file in ([options.configfile]):
		if not (file):
			parser.print_help()
			sys.exit(0)

	config_json = json.loads(open(options.configfile, "r").read())

	progress_file = config_json["mutmapper"]["outdir"] + "/mutmap-step2-progress.txt"
 	mutmap_file = config_json["mutmapper"]["outdir"] + "/mutmap.1.csv"
	out_file = config_json["mutmapper"]["outdir"] + "/mutmap.2.csv"



	is_canon = {}
	with open(config_json["idmapping"], "r") as FR:
                data_grid = csv.reader(FR, delimiter=',', quotechar='"')
                row_count = 0
                for row in data_grid:
                        isoform = row[0]
        		is_canon[isoform] = True 


	pepid2acmap = {"uniprot1":{}, "refseq1":{}, "uniprot2":{}, "refseq2":{}}
	with open(config_json["pepid2isoformac"], "r") as FR:
                data_grid = csv.reader(FR, delimiter=',', quotechar='"')
                row_count = 0
                for row in data_grid:
                        pep_id, ac_list = row[0], row[2].split("|")
                        if row[1] == "exact":
				pepid2acmap["uniprot1"][pep_id] = "|".join(ac_list)

	with open(config_json["pepid2refseqac"], "r") as FR:
                data_grid = csv.reader(FR, delimiter=',', quotechar='"')
                row_count = 0
                for row in data_grid:
                        pep_id, ac_list = row[0], row[2].split("|")
                        if row[1] == "exact":
                                pepid2acmap["refseq1"][pep_id] = "|".join(ac_list)


	idmap_file = config_json["mutmapper"]["outdir"] + "/map2uniprot.ranges.csv"
	with open(idmap_file, "r") as FR:
                data_frame = csv.reader(FR, delimiter=',', quotechar='"')
                row_count = 0
		for row in data_frame:
                        row_count += 1
			if row_count == 1:
				continue
			pep_id, ac = row[0], row[1]
			pepid2acmap["uniprot2"][pep_id] = {}
			for r in row[2].strip().split("|"):
				r1, r2 = r.split(":")
				start_pos1,end_pos1 = r1.split("-")
				start_pos2,end_pos2 = r2.split("-")
				n = 0
				for pos in xrange(int(start_pos1), int(end_pos1) + 1):
					val = "%s %s" % (ac,int(start_pos2) + n)
					pepid2acmap["uniprot2"][pep_id][str(pos)] = val
					n += 1

	idmap_file = config_json["mutmapper"]["outdir"] + "/map2refseq.ranges.csv"
        with open(idmap_file, "r") as FR:
                data_frame = csv.reader(FR, delimiter=',', quotechar='"')
		row_count = 0
                for row in data_frame:
                        row_count += 1
                        if row_count == 1:
                                continue
                        pep_id, ac = row[0], row[1]
                        pepid2acmap["refseq2"][pep_id] = {}
                        for r in row[2].strip().split("|"):
				r1, r2 = r.split(":")
                                start_pos1,end_pos1 = r1.split("-")
                                start_pos2,end_pos2 = r2.split("-")
                                n = 0
                                for pos in xrange(int(start_pos1), int(end_pos1) + 1):
                                        val = "%s %s" % (ac,int(start_pos2) + n)
					pepid2acmap["refseq2"][pep_id][str(pos)] = val
                                        n += 1

	
	FW = open(out_file, "w")
	with open(mutmap_file, 'rb') as FR:
		data_frame = csv.reader(FR, delimiter=',', quotechar='"')
		row_count = 0
		for row in data_frame:
                        row_count += 1
			if row_count == 1:
				continue
			pep_id = row[6]
			position_in_pep = row[8]
			uniprot_ac, pos_in_uniprot, refseq_ac, pos_in_refseq,flag1,flag2 = "-", "-", "-", "-", "-","-"
			if pep_id in pepid2acmap["uniprot2"]:
				if position_in_pep in pepid2acmap["uniprot2"][pep_id]:
					uniprot_ac = pepid2acmap["uniprot2"][pep_id][position_in_pep].split(" ")[0]
					pos_in_uniprot = pepid2acmap["uniprot2"][pep_id][position_in_pep].split(" ")[1]
					flag1 = "nonexact"
			if pep_id in pepid2acmap["refseq2"]:
                                if position_in_pep in pepid2acmap["refseq2"][pep_id]:
                                        refseq_ac = pepid2acmap["refseq2"][pep_id][position_in_pep].split(" ")[0]
                                        pos_in_refseq = pepid2acmap["refseq2"][pep_id][position_in_pep].split(" ")[1]
					flag2 = "nonexact"
			if pep_id in pepid2acmap["uniprot1"]:
				uniprot_ac = pepid2acmap["uniprot1"][pep_id]
				pos_in_uniprot = position_in_pep
				flag1 = "exact"
			if pep_id in pepid2acmap["refseq1"]:
                                refseq_ac = pepid2acmap["refseq1"][pep_id]
                                pos_in_refseq = position_in_pep
 				flag2 = "exact"
			newrow = row + [uniprot_ac, pos_in_uniprot,flag1, refseq_ac, pos_in_refseq,flag2,
					str(uniprot_ac in is_canon)]
			FW.write("%s\n" % (",".join(newrow)))
	FW.close()

if __name__ == '__main__':
        main()








