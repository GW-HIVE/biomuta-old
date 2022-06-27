import os,sys
import string
import csv
import json
from optparse import Option_parser
from Bio import Seq_iO
from Bio.Seq import Seq


__version__="1.0"
__status__ = "Dev"


#################################
def report_progress(progress_file, msg, open_new=False):

        if open_new == True:
                with open(progress_file, "w") as FW:
                        FW.write(msg)
        else:
                with open(progress_file, "a") as FA:
                        FA.write(msg)
        return


###################
def get_attr_dict(attr_string):

        att_hash = {}
        attr_list = attr_string.strip().split(";")
        for attr in attr_list:
                if attr != "":
                        name, value = attr.strip().split(" ")
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
def load_transcript_objects(gtf_file, chr_id_list, genome_seq_hash, progress_file):


        
        
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

                                	if row_count%25000 == 0:
						msg = "Loaded %s exons from GTF file\n" % (row_count)
                                        	report_progress(progress_file, msg)
                                	row_count += 1


	return exon_objs




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
	chr_list = map(str,range(1,23)) + ["X","Y","MT"]
	#chr_list = ["1"]

	cds_file = config_json["cdsfile"]
	gtf_file = config_json["gtffile"]
	genome_file = config_json["genomefile"]
	exon_list_file = config_json["mutmapper"]["outdir"] + "/exonlist.json"
	transcript_list_file = config_json["mutmapper"]["outdir"] + "/transcriptlist.csv"
	progress_file = config_json["mutmapper"]["outdir"] + "/cdscheck-progress.txt"
        


	
	#Dump transcript list
        FW = open(transcript_list_file, "w")
	report_progress(progress_file, "Started dumping transcript list\n", True)
	with open(gtf_file, 'rb') as tsvfile:
                tsvreader = csv.reader(tsvfile, delimiter='\t', quotechar='|')
                row_count = 0
                for row in tsvreader:
                        if len(row) > 2:
                                c1 = row[2] in ["transcript"]
                                if c1:
					att_hash = get_attr_dict(row[-1])
                        		FW.write("%s,%s\n" % (att_hash["transcript_id"], row[1]))
	FW.close()
 



	#Load genomic sequences
	genome_seq_hash = {}
        for record in Seq_iO.parse(genome_file, "fasta"):
                chr_id = record.id.strip().replace("chr", "")
		if chr_id in chr_list:
			genome_seq_hash[chr_id] = record.seq.upper()
			with open(progress_file, "a") as FA:
                		msg = "loaded genomic sequences for chr = %s\n" %(chr_id)
				FA.write(msg)
               

	#Load cds sequences
	report_progress(progress_file, "\n_started loading CDS sequences ...\n" )
	cds_seq_hash = {}
	for record in Seq_iO.parse(cds_file, "fasta"):
		cds_seq_hash[record.id] = record.seq.upper()
	report_progress(progress_file, "Done.\n" )
        

	#Load transcript objects
	report_progress(progress_file, "\n_started loading exon objects from GTF file ...\n" )
	exon_objs = load_transcript_objects(gtf_file, chr_list, genome_seq_hash, progress_file) 
	report_progress(progress_file, "Done.\n" )


	flags_file = config_json["mutmapper"]["outdir"] + "/cds_ranges_flags.txt"
       	FW1 = open(flags_file, "w")
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
			if seq1 == seq2:
				FW1.write("1|%s|%s\n" % ("1", trs_id))
			if seq1 != seq2:
				flag, seq2 = correct_harvested_cds(seq2,frame1, frame2)
				if seq1 == seq2:
					exon_objs[trs_id]["correctionflag"] = flag
					FW1.write("1|%s|%s\n" % (flag,trs_id))
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
						
				else:
					FW1.write("0|%s|%s\n" % ("x",trs_id))	
		else:
			FW1.write("0|%s|%s\n" % ("y",trs_id))

	FW1.close()
	


	report_progress(progress_file, "\n_now writing out exon objects to JSON file ...\n")
	FW = open(exon_list_file, "w")
	FW.write("%s\n" % (json.dumps(exon_objs, indent=4)))
	FW.close()
	report_progress(progress_file, "Done.\n")


	report_progress(progress_file, "\n_now loading exon objects from file and checking exon to cds coordinate mapping ... \n")
	exon_objs = json.loads(open(exon_list_file, "r").read())
        for trs_id in exon_objs:
                flag = True
                for obj in exon_objs[trs_id]["exonlist"]:
                        if obj["status"] == "removed":
                                continue
                        cds_start, cds_stop = obj["cdsstart"] - 1, obj["cdsstop"] - 1
                        seq1 = cds_seq_hash[trs_id][cds_start:cds_stop]
                        seq2 = obj["seq"]
                        flag = False if seq1 != seq2 else flag
                if flag == False:
			msg = "Check failed for %s %s %s\n" %(trs_id, flag, exon_objs[trs_id]["correctionflag"])
                        report_progress(progress_file, msg)
                        sys.exit()
	report_progress(progress_file, "All exons, except those removed, map to cds sequences correctly\n")







if __name__ == '__main__':
        main()








