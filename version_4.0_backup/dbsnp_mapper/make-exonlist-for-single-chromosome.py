import os,sys
import string
import csv
import json
from Bio import SeqIO
from Bio.Seq import Seq
from optparse import OptionParser



__version__="1.0"
__status__ = "Dev"




###################
def get_attr_dict(attr_string):


    att_hash = {}
    attr_list = attr_string.strip().split(";")
    for attr in attr_list:
        if attr != "":
            name = attr.strip().split(" ")[0]
            value = " ".join(attr.strip().split(" ")[1:])
            att_hash[name] = value.replace("\"", "")

    return att_hash


###########################
def correct_harvested_cds(harvested_seq, frame1, frame2):

	if frame1 == 2:
		return "a", "N" + harvested_seq
	elif frame1 == 1:
		return "b", "NN" + harvested_seq
	elif frame2 == 1:
                return "c", harvested_seq[:-4] + harvested_seq[-3:]
	elif frame2 == 2:
                return "d", harvested_seq[:-5] + harvested_seq[-3:]
	else:
		return "x", harvested_seq






########################################
def load_transcript_objects(gtf_file, chr_id_list, genome_seq_hash):


        
        
	exon_objs = {}
        with open(gtf_file, 'rb') as tsvfile:
                tsvreader = csv.reader(tsvfile, delimiter='\t', quotechar='|')
                row_count = 0
                for row in tsvreader:
                        if len(row) > 2:
				cond_list = []
				cond_list.append(row[0] in chr_id_list)
				cond_list.append(row[2] in ["CDS", "stop_codon"])
				if False not in cond_list:
					att_hash = get_attr_dict(row[-1])
					obj = {
						"cds_type":row[2],
						"chr_id":row[0],
						"gene_name":att_hash["gene_name"],
						"transcript_id":att_hash["transcript_id"],
						"frame":int(row[7]),
						"exon_number":int(att_hash["exon_number"]),
						"start":int(row[3]), "stop":int(row[4]),
						"strand":row[6]
					}
					chr_id = obj["chr_id"]
					trs_id = obj["transcript_id"]
					i, j = int(obj["start"]) - 1, int(obj["stop"])
                                        obj["seq"] = ""
                                        if chr_id in genome_seq_hash:
                                                obj["seq"] = genome_seq_hash[chr_id][i:j]
                                                if obj["strand"].strip() == "-":
                                                        tmp_seq = Seq(str(obj["seq"]))
                                                        obj["seq"] =  tmp_seq.reverse_complement()
                                        obj["seq"] = str(obj["seq"])
					if trs_id not in exon_objs:
                                                exon_objs[trs_id] = {"exonlist":[]}
                                        exon_objs[trs_id]["exonlist"].append(obj)

                                	row_count += 1


	return exon_objs




###############################
def main():



        gtf_file = "/data/external/ucsc/grch37/gtf/Homo_sapiens.GRCh37.75.gtf"
        dna_seq_file = "/data/external/ucsc/grch37/targets/Homo_sapiens.GRCh37.75.dna.all.fa"
        cds_seq_file = "/data/external/ucsc/grch37/targets/Homo_sapiens.GRCh37.75.cds.all.fa"


	chr_list = ["12"]

	#Load genomic sequences
	genome_seq_hash = {}
        for record in SeqIO.parse(dna_seq_file, "fasta"):
                chr_id = record.id.strip().replace("chr", "")
		if chr_id in chr_list:
			genome_seq_hash[chr_id] = record.seq.upper()
               

	#Load cds sequences
	cds_seq_hash = {}
	for record in SeqIO.parse(cds_seq_file, "fasta"):
	    seq_id = record.id.split(".")[0]
            cds_seq_hash[seq_id] = record.seq.upper()
        

	#Load transcript objects
	exon_objs = load_transcript_objects(gtf_file, chr_list, genome_seq_hash) 

	for trs_id in exon_objs:
		exon_objs[trs_id]["correctionflag"] = "none"
		frame1 = exon_objs[trs_id]["exonlist"][0]["frame"]
		frame2 = exon_objs[trs_id]["exonlist"][-1]["frame"]
		cds_type = exon_objs[trs_id]["exonlist"][-1]["cds_type"]
		seq2 =  ""
		cds_len = 1
		for obj in exon_objs[trs_id]["exonlist"]:
			seq2 += obj["seq"]
			exon_len = obj["stop"] - obj["start"] + 1
			obj["cdsstart"], obj["cdsstop"], obj["status"] = cds_len, cds_len + exon_len, ""
			cds_len += exon_len
		if trs_id in cds_seq_hash:
			seq1 = cds_seq_hash[trs_id]
			if seq1 != seq2:
				flag, seq2 = correct_harvested_cds(seq2,frame1, frame2)
				if seq1 == seq2:
					exon_objs[trs_id]["correctionflag"] = flag
					if flag in ["a", "b"]:
						incr = 1 if flag == "a" else 2
						for obj in exon_objs[trs_id]["exonlist"]:
							obj["cdsstart"] += incr
							obj["cdsstop"] += incr
					elif flag in ["c", "d"]:
						incr = 1 if flag == "c" else 2
						exon_objs[trs_id]["exonlist"][-2]["cdsstart"] -= incr
						exon_objs[trs_id]["exonlist"][-2]["cdsstop"] -= incr
						exon_objs[trs_id]["exonlist"][-1]["cdsstart"] -= incr
                                                exon_objs[trs_id]["exonlist"][-1]["cdsstop"] -= incr
						exon_objs[trs_id]["exonlist"][-3]["status"] = "removed"
	
	print json.dumps(exon_objs, indent=4)









if __name__ == '__main__':
        main()








