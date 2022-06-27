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

	genome_file = config_json["genomefile"]
	cds_file = config_json["cdsfile"]
	pep_file = config_json["pepfile"]
	codon_file = config_json["codonfile"]
	
	exon_list_file = config_json["mutmapper"]["outdir"] + "/exonlist.json"
	mutlist_file = config_json["mutmapper"]["outdir"] + "/mutlist-merged.csv"
	progress_file = config_json["mutmapper"]["outdir"] + "/mutmap-step1-progress.txt"
 	mutmap_file = config_json["mutmapper"]["outdir"] + "/mutmap.1.csv"
	flags_file = config_json["mutmapper"]["outdir"] + "/map_flags.txt"
	map_stat_file = config_json["mutmapper"]["outdir"] + "/map_stat.csv"


	chr_list = map(str,range(1,23)) + ["X","Y","MT"]
     
	#chr_list = ["17"]


	#Load genomic sequences
	report_progress(progress_file, "\n_started loading genome sequences ...\n", True)
        genome_seq_hash = {}
        for record in Seq_iO.parse(genome_file, "fasta"):
                chr_id = record.id.strip().replace("chr", "")
                if chr_id in chr_list:
                        genome_seq_hash[chr_id] = record.seq.upper()
	report_progress(progress_file, "Done.\n")



	  
	#Load cds sequences
        report_progress(progress_file, "\n_started loading CDS sequences ...\n")
        cds_seq_hash = {}
        for record in Seq_iO.parse(cds_file, "fasta"):
                cds_seq_hash[record.id] = record.seq.upper()
        report_progress(progress_file, "Done.\n")


	#Load pep sequences
        report_progress(progress_file, "\n_started loading peptide sequences ...\n")
        pep_seq_hash = {}
	trsid2pepid = {}
	i = 0
        for record in Seq_iO.parse(pep_file, "fasta"):
                desc = record.description
                trs_id = desc.split("transcript:")[1].split()[0]
                pep_seq_hash[trs_id] = str(record.seq.upper())
		trsid2pepid[trs_id] =  record.id
                i += 1
        report_progress(progress_file, "Done.\n")




	#Load codon dict
        report_progress(progress_file, "\n_started loading codon dictionary ...\n")
        codon_dic = {}
        with open(codon_file, "r") as tsvfile:
                for row in tsvfile:
                        row = row.strip().split('\t')
                        codon_dic[row[0]] = row[1]
	report_progress(progress_file, "Done.\n")


	report_progress(progress_file, "\n_started loading exon objects from JSON file ...\n")
	exon_objs = json.loads(open(exon_list_file, "r").read())
	report_progress(progress_file, "Done.\n")


	#Populate pos to trsid mapping
	report_progress(progress_file, "\n_started constructing cdslocus2tridlist map ...\n")
	cdslocus2tridlist = {}
	for trs_id in exon_objs:
                for obj in exon_objs[trs_id]["exonlist"]:
                        for j in xrange(int(obj["start"]), int((obj["stop"])) + 1):
                                k = obj["chr_id"] + ":" + str(j)
				if k not in cdslocus2tridlist:
					cdslocus2tridlist[k] = []
                                if trs_id not in cdslocus2tridlist[k]:
                                        cdslocus2tridlist[k].append(trs_id)
	report_progress(progress_file, "Done.\n")




	report_progress(progress_file, "\n_started loading mutation list ...\n")
        mut_list, mut_count_hash  = load_mut_list_in_cds(mutlist_file, chr_list, cdslocus2tridlist)
        #mut_list = [[17,41246211,"C","G","True","TCGA","-"]]
	report_progress(progress_file, "Done.\n")

	FW = open(mutmap_file, "w")
	report_progress(flags_file, "Flags", True)


	cmp_hash = {"A":"T", "T":"A", "C":"G", "G":"C"}
	report_progress(progress_file, "\n_started identifying mutation effect on translation ...\n")
        for mut_row in mut_list:
		mut_index,chr_id, mut_pos_in_genome, ref, alt = mut_row[0],mut_row[1],int(mut_row[2]),mut_row[3],mut_row[4]
		debug_flag = chr_id == "7" and mut_pos_in_genome == 142448658
		locus = mut_row[1] + ":" + mut_row[2] 
		for trs_id in cdslocus2tridlist[locus]:
			if debug_flag:
				print "flag-1", trs_id
			frame1 = exon_objs[trs_id]["exonlist"][0]["frame"]
			frame2 = exon_objs[trs_id]["exonlist"][-1]["frame"]
                        for exon_index in xrange(0, len(exon_objs[trs_id]["exonlist"])):
                                obj = exon_objs[trs_id]["exonlist"][exon_index]
				alt_new = cmp_hash[alt] if obj["strand"] == "-" else alt
				ref_new = cmp_hash[ref] if obj["strand"] == "-" else ref
				#disregard removed exons
				if obj["status"] == "removed":
					if debug_flag:
						print "flag-x"
					continue
				if mut_pos_in_genome >= int(obj["start"]) and mut_pos_in_genome <=  int(obj["stop"]):
					if debug_flag:
		                                print "flag-2", obj["start"], obj["stop"]

					exon_len = obj["stop"] - obj["start"] + 1
                                        mut_pos_in_exon = mut_pos_in_genome - obj["start"]
                                        if obj["strand"] == "-":
                                                mut_pos_in_exon = obj["stop"] - mut_pos_in_genome
                                        mut_pos_in_cds = obj["cdsstart"] - 1 + mut_pos_in_exon
                                        mut_pos_in_codon = (mut_pos_in_exon-obj["frame"])%3
                                        ref_in_genome = genome_seq_hash[chr_id][mut_pos_in_genome-1]
					ref_in_exon = obj["seq"][mut_pos_in_exon]
                                        ref_in_cds = cds_seq_hash[trs_id][mut_pos_in_cds]
					codon_index_in_cds = int(mut_pos_in_cds/3)
					ind = codon_index_in_cds*3
					codon_seq = str(cds_seq_hash[trs_id][ind:ind+3])


					last_codon_flag = codon_index_in_cds == len(pep_seq_hash[trs_id]) 
					if debug_flag:
                                                print "flag-3"
					aa_translated = translate_codon(codon_seq, codon_dic, last_codon_flag)
					if debug_flag:
                                                print "flag-4"
                                        	
					out_values = []
					flag = ""
					if ref_in_exon != ref_in_genome and ref_in_exon != cmp_hash[ref_in_genome]:
						flag = "G"
						out_values = [flag]+ mut_row[0:5] 
						out_values += [trs_id,str(mut_pos_in_genome),ref_in_exon,ref_in_genome]
					elif ref_in_exon != codon_seq[mut_pos_in_cds%3]:
                                                flag = "A"
                                                out_values = [flag]+ mut_row[0:5]
						out_values += [trs_id,ref_in_exon,codon_seq,str(mut_pos_in_cds%3)]
                                        elif ref_in_exon != ref_in_cds:
                                                flag = "B"
                                                out_values = [flag]+ mut_row[0:5]
						out_values += [trs_id,str(mut_pos_in_genome),str(mut_pos_in_exon),
								ref_in_exon,
								ref_in_cds,obj["status"]]
					elif aa_translated != "":
						if debug_flag:
                                        	        print "flag-5"
                                        
						if codon_index_in_cds < len(pep_seq_hash[trs_id]):
							if debug_flag:
                                        		        print "flag-6"
							aa_in_pep_seq = pep_seq_hash[trs_id][codon_index_in_cds]	
							if aa_translated == aa_in_pep_seq:
								if debug_flag:
                                        			        print "flag-7"
								alt_codon_seq = codon_seq[:mut_pos_in_codon] + alt_new + codon_seq[mut_pos_in_codon+1:]
								aa_alt = translate_codon(alt_codon_seq,
									codon_dic, False)
								if alt_codon_seq == codon_seq:
									if debug_flag:
                                        				        print "flag-8"
									flag = "fatal-H"
									out_values = [flag,str(mut_index),
									trs_id,
									ref_new,alt_new,codon_seq,
									alt_codon_seq,str(mut_pos_in_codon)]
								else:
									if debug_flag:
                                                				print "flag-9"
									out_row = mut_row[0:5]
									out_row.append(trs_id)
									out_row.append(trsid2pepid[trs_id])
									out_row.append(str(mut_pos_in_cds + 1))
									out_row.append(str(codon_index_in_cds + 1))
									out_row.append(str(mut_pos_in_codon + 1))
									out_row.append(codon_seq)
									out_row.append(alt_codon_seq)
									out_row.append(aa_translated)
									out_row.append(aa_alt)
									line_out = ",".join(out_row)	
									FW.write("%s\n" % (line_out))
							else:
								if debug_flag:
                                       				         print "flag-10"
								flag = "C"
								out_values = [flag] + mut_row[0:5]
								out_values += [trs_id,codon_seq,
										aa_translated,aa_in_pep_seq]
						else:
							if debug_flag:
                                                		print "flag-11"
                                        		flag = "D"
							out_values = [flag] + mut_row[0:5]
							out_values += [trs_id,codon_seq,
									aa_translated,
									codon_index_in_cds, 
									len(pep_seq_hash[trs_id])]
					else:
						if debug_flag:
                                                	print "flag-12"
                                        	flag = "E"
						out_values = [flag] + mut_row[0:5]
						out_values += [trs_id,codon_seq,
								aa_translated,
								aa_in_pep_seq,str(last_codon_flag)]
					if flag != "":
						if debug_flag:
                                        	        print "flag-13"
						for ii in xrange(0, len(out_values)):
							out_values[ii] = str(out_values[ii])
                                                msg = ",".join(out_values)
                                                report_progress(flags_file, msg)

	FW.close()

	
	with open(map_stat_file, "w") as FW:
		for key in ["total", "coding", "noncoding", "incontigs"]:
			FW.write("%s,%s\n" %(key, mut_count_hash[key]))
	

	report_progress(progress_file, "Successfully finished.\n")




if __name__ == '__main__':
        main()








