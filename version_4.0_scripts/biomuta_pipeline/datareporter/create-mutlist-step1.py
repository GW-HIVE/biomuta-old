import os,sys
import My_sQLdb
import string
import commands
from optparse import Option_parser
import glob
import json
import csv


__status__ = "Dev"
__version__ = "4.0"


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




##########################################
def parse_line_single(src_type, vcf_file, row_count, values):


	pass_filter = True
	combo_id_list = []
        extra_cols = []
	chr_id = values[0] #Get single nucleotide mutation

	if src_type == "tcga":
		if values[6].strip() != "PASS":
                        pass_filter = False
		cancer_type = vcf_file.split("/")[-3]
		combo_id_list.append(values[0] + "," + values[1] + "," + values[3] + "," + values[4])
		extra_cols.append(str(pass_filter)+ ",TCGA")
	elif src_type == "cosmic":
		if "SNP" in values[-1].split(";"):
                        pass_filter = False
		if values[2] in cosmicid2doid:
			cancer_type = cosmicid2doid[values[2]]
		else:
			cancer_type = values[2]
			pass_filter = False
		combo_id_list.append(values[0] + "," + values[1] + "," + values[3] + "," + values[4])
		if vcf_file.find("WGS") >= 0:
			extra_cols.append(str(pass_filter) + ",COSMIC-WGS")
		elif vcf_file.find("Coding") >= 0:
			extra_cols.append(str(pass_filter) + ",COSMIC-WXS")
	elif src_type == "clinvar":
		if values[-1].split(';')[4] not in ["SAO=0", "SAO=2", "SAO=3"]:
                        pass_filter = False
		cancer_type = values[-1].split(";CLNDBN=")[-1].split(';')[0]
                combo_id_list.append(values[0]+","+values[1] + "," + values[3] + "," + values[4])
		extra_cols.append(str(pass_filter) + ",CLINVAR")	
	elif src_type == "icgc":
		cancer_type = values[-1].split(";OCCURRENCE=")[-1].split("-")[0]
		combo_id_list.append(values[0] + "," + values[1] + "," + values[3] + "," + values[4])
		extra_cols.append(str(pass_filter) + ",ICGC")
	elif src_type == "civic":
		cancer_type = values[-2]
		combo_id_list.append(values[0] + "," + values[1] + "," + values[2] + "," + values[3])
                extra_cols.append(str(pass_filter) + ",CIVIC")	
	
	return cancer_type, combo_id_list, extra_cols



#######################################
def load_do_slim (doid_file, doslimid_hash, doslimname_hash, exclude_list):


        with open(doid_file, 'r') as FR:
                csv_grid = csv.reader(FR, delimiter='\t', quotechar='"')
                for row in csv_grid:
                        doslim_id = row[0].split(" ")[0].split(":")[1]
                        if doslim_id in exclude_list:
                                continue
                        doslimname_hash[doslim_id] = row[0]
                        doslimid_hash[row[0]] = doslim_id
                        doslimid_hash[doslim_id] = doslim_id
			for val in row[1].split("|"):
                                doslimid_hash[val] = doslim_id
                        for j in [2,3]:
                                val = row[j].strip()
                                if val != "":
                                        doslimid_hash[val] = doslim_id
        return



##############################################
def main():

        usage = "\n%prog  [options]"
        parser = Option_parser(usage,version="%prog " + __version__)
        parser.add_option("-i", "--configfile", action = "store", dest = "configfile", help = "Config json file path")
	parser.add_option("-s","--srctype",action="store",dest="srctype",help="Source type")
        
	(options,args) = parser.parse_args()
        for file in ([options.srctype]):
                if not (file):
                        parser.print_help()
                        sys.exit(0)

	src_type = options.srctype
	json_file = options.configfile
        config_json = json.loads(open(json_file).read())

	input_dir = config_json["datareporter"][src_type]["input_dir"]
	output_dir = config_json["datareporter"][src_type]["output_dir"]

	global count_hash
	global seen_hash
	global cosmicid2doid

		
	count_hash = {
		"badrows":{},
		"variants_badrows":0,
		"variants_total":0,
		"variants_valid":0,
		"variants_valid_unique":0
	}
	doslimid_hash = {}
	doslimname_hash = {}
	cosmicid2doid = {}
	
	cosmicdoid_file = output_dir + "/cosmic-mutid2doid.tsv"
	if src_type == "cosmic":
		with open(cosmicdoid_file, 'r') as FR:       
 			csv_grid = csv.reader(FR, delimiter='\t', quotechar='"')
			for row in csv_grid:
				cosmicid2doid[row[0]] = row[1]



	#Load bad rows
	bad_row_hash = {}	
	bad_rows_file = output_dir + "/badrows.txt"
	with open(bad_rows_file, 'r') as FR:
		csv_grid = csv.reader(FR, delimiter=',', quotechar='"')
		for row in csv_grid:
			if row[0] not in count_hash["badrows"]:
				count_hash["badrows"][row[0]] = 0	
			count_hash["badrows"][row[0]] += 1
			if row[2] not in bad_row_hash:
				bad_row_hash[row[2]] = {}
			bad_row_hash[row[2]][row[1]] = True

	#Load DOID mappings
	doid_file = config_json["datareporter"]["doid"]["output_dir"] + "/doid-mapping.reduced.tsv"
	load_do_slim (doid_file, doslimid_hash, doslimname_hash, [])
        #Repeat loading to overwrite mapping to 305 if more specific mapping exists
	load_do_slim (doid_file, doslimid_hash, doslimname_hash, ["305"])




	chr_list1 = map(str,range(1,23)) + ["X","Y"]
	chr_list2 = ["1"]
	if src_type == "tcga":
		pattern = input_dir + "/*/*/*.vcf"
	else:
		pattern = input_dir + "/*.vcf"
	file_list = glob.glob(pattern)

	#Initialize progress file	
	progress_file = output_dir + "/create-mutlist-progress.txt"
        with open(progress_file, "w") as FW:
        	FW.write("Started writing to progress file\n")

	
	nx = 0
	cancer_type_hash = {}
	chr_list = chr_list2 if src_type == "tcga" else chr_list1
	for chr_id in chr_list:
		freq_hash = {}
        	seen_hash = {"comboid":{}, "nondoslim":{}}
		newdoslim_id = 0
		file_count = 0
		for vcf_file in file_list:
			cancer_type = vcf_file.split("/")[-3]
			file_name = vcf_file.split("/")[-1]
			file_count += 1
			with open(vcf_file, 'r') as FR:
                		row_count = 0
				for line in FR:
					row_count += 1
					row = line.strip().split("\t")
					if row[0][0:1] == "#":	
						continue
					if src_type!= "tcga" and row[0] != chr_id:
						continue
					chr_id = row[0] if src_type == "tcga" else chr_id 
					if chr_id not in chr_list1:
						continue                                        
					count_hash["variants_total"] += 1
					if vcf_file in bad_row_hash:
						if str(row_count) in bad_row_hash[vcf_file]:
							count_hash["variants_badrows"] += 1
							continue 
					cancer_type, combo_id_list, extra_cols_list = parse_line_single(src_type,vcf_file, 
											row_count, row)
					
					cancer_type_hash[cancer_type] = True
					if cancer_type not in doslimid_hash:
						newdoslim_id -= 1
						doslimid_hash[cancer_type] = str(newdoslim_id)
						doslimname_hash[str(newdoslim_id)] = "NONDOSLIM"
						seen_hash["nondoslim"][cancer_type] = True
					doslim_id = doslimid_hash[cancer_type]	

					if len(combo_id_list) > 0:
						for j in xrange(0, len(combo_id_list)):
							combo_id = combo_id_list[j]
							extra_cols = extra_cols_list[j]
							if combo_id not in seen_hash["comboid"]:
								seen_hash["comboid"][combo_id] = []
								count_hash["variants_valid_unique"] += 1
							long_combo_id = combo_id + "," + extra_cols
							if long_combo_id not in freq_hash:
								freq_hash[long_combo_id] = {}
							if doslim_id not in freq_hash[long_combo_id]:
								freq_hash[long_combo_id][doslim_id] = 1
							else:
								freq_hash[long_combo_id][doslim_id] += 1
							count_hash["variants_valid"] += 1
				
					if row_count%1000000 == 0:
						report_progress(progress_file,"chr%s: parsed %s rows" % (chr_id,row_count))
			report_progress(progress_file, "processed %s/%s files.\n" % (file_count, len(file_list)))
	
		chr_lbl = chr_id if src_type != "tcga" else "all"
		mut_list_file = output_dir + "/mutlist.chr"+chr_lbl+".csv"
		with open(mut_list_file, "w") as FW:
                	FW.write("chr,pos,ref,alt,passfilter,source,frequency\n")
               	 	for long_combo_id in freq_hash:
                        	values = []
                        	for doslim_id in freq_hash[long_combo_id]:
                               		string = doslim_id + ":" + str(freq_hash[long_combo_id][doslim_id])
					values.append(string)
                        	FW.write("%s,%s\n" % (long_combo_id," ".join(values)))





	nondoslim_file = output_dir + "/doslim.tsv"
	with open(nondoslim_file, "w") as FW:
		FW.write("doslimid\tdoslimname\tcancertype\n")
		#for cancer_type in seen_hash["nondoslim"]:
		for cancer_type in cancer_type_hash:
			doslim_id = doslimid_hash[cancer_type]
			doslim_name = doslimname_hash[doslim_id]
			FW.write("%s\t%s\t%s\n" % (doslim_id,doslim_name,cancer_type))

	stat_file = output_dir + "/vcf.stat.csv"
	with open(stat_file, "w") as FW:
		FW.write("count,type\n")
		FW.write("%s,%s\n" % (count_hash["variants_total"], "total number of calls"))
		FW.write("%s,%s\n" % (count_hash["variants_badrows"], "problematic calls (badrows)"))
		FW.write("%s,%s\n" % (count_hash["variants_valid"], "valid calls"))
		FW.write("%s,%s\n" % (count_hash["variants_valid_unique"], "valid unique calls"))
			

	stat_file = output_dir + "/badrows.stat.csv"
        with open(stat_file, "w") as FW:
		FW.write("count,type\n")
		for k in sorted(count_hash["badrows"], key=lambda x: count_hash["badrows"][x], reverse=True):
                       	FW.write("%s,%s\n"%(count_hash["badrows"][k],k))
		FW.write("%s,%s\n"%(count_hash["variants_badrows"], "total"))






if __name__ == '__main__':
	main()


