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

	msg = msg.strip()
        if open_new == True:
                with open(progress_file, "w") as FW:
                        FW.write(msg + "\n")
        else:
                with open(progress_file, "a") as FA:
                        FA.write(msg + "\n")
        return

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
                        for val in row[1].split("|"):
                                doslimid_hash[val] = doslim_id
                        for j in [2,3]:
                                val = row[j].strip()
                                if val != "":
                                        doslimid_hash[val] = doslim_id
        return




###############################
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

	
	chr_list = map(str,range(1,23)) + ["X","Y"]
	
	count_hash = {"pass":{}, "all":{}}
	for chr_id in chr_list:
		mutlist_file = config_json["datareporter"][src_type]["output_dir"]
		mutlist_file += "/mutlist.chr"+chr_id+".csv"
		if os.path.exists(mutlist_file) == False:
			continue
		with open(mutlist_file, 'rb') as tsvfile:
                       	tsvreader = csv.reader(tsvfile, delimiter=',', quotechar='"')
                       	row_count = 0
			for row in tsvreader:
                               	row_count += 1
				if row_count == 1:
                                       	continue
				if row[0] != chr_id:
					continue 
				pass_filter = row[4]
				combo_id = ",".join([row[0], row[1], row[2], row[3]])
                               	term_list = row[6].strip().split(" ")
				for term in term_list:
					doslim_id = term.split(":")[0]
					freq = int(term.split(":")[1])
					if doslim_id not in count_hash["all"]:
						count_hash["all"][doslim_id] = {}
					if chr_id not in count_hash["all"][doslim_id]:
						count_hash["all"][doslim_id][chr_id] = freq
					else:
						count_hash["all"][doslim_id][chr_id] += freq
					if pass_filter == "True":
						if doslim_id not in count_hash["pass"]:
                                                	count_hash["pass"][doslim_id] = {}
                                        	if chr_id not in count_hash["pass"][doslim_id]:
                                                	count_hash["pass"][doslim_id][chr_id] = freq
                                        	else:
                                                	count_hash["pass"][doslim_id][chr_id] += freq	


	doslimid_hash = {}
        doslimname_hash = {"total":"total"}
	cancertype_hash = {}
	doslimmap_file = config_json["datareporter"][src_type]["output_dir"] + "/doslim.tsv"
	with open(doslimmap_file, 'rb') as tsvfile:
		tsvreader = csv.reader(tsvfile, delimiter='\t', quotechar='"')
		row_count = 0
		for row in tsvreader:
			row_count += 1
			if row_count == 1:
				continue
			doslim_id, doslim_name, cancer_type = row[0], row[1], row[2]
			if int(doslim_id) < 0:
				cancertype_hash[doslim_id] = cancer_type
			else:
				doslimname_hash[doslim_id] = doslim_name

	
	stat_file = config_json["datareporter"][src_type]["output_dir"] + "/allmuts.stat.csv"
        FW = open(stat_file, "w")

        FW.write("%s\n" % ("doslimname,chr" + ",chr".join(chr_list) + ",total"))
	doslim_list = count_hash["all"].keys() + ["total"] 
	count_hash["all"]["total"] = {}
	for doslim_id in doslim_list:
		count_hash["all"][doslim_id]["total"] = 0
		doslim_name = doslimname_hash[doslim_id] if doslim_id in doslimname_hash else cancertype_hash[doslim_id] 
		
		values = ["\""+doslim_name+"\""]
		for chr_id in chr_list + ["total"]:
			val = count_hash["all"][doslim_id][chr_id] if chr_id in count_hash["all"][doslim_id] else 0
			values.append(str(val))
			count_hash["all"][doslim_id]["total"] += val
			if chr_id in count_hash["all"]["total"]:
				count_hash["all"]["total"][chr_id] += val
			else:
				count_hash["all"]["total"][chr_id] = val
		FW.write("%s\n" % (",".join(values)))
	FW.close()


	stat_file = config_json["datareporter"][src_type]["output_dir"] + "/passedmuts.stat.csv"
       	FW = open(stat_file, "w")
        FW.write("%s\n" % ("doslimname,chr" + ",chr".join(chr_list) + ",total"))
        doslim_list = count_hash["pass"].keys() + ["total"]
        count_hash["pass"]["total"] = {}
        for doslim_id in doslim_list:
                count_hash["pass"][doslim_id]["total"] = 0
                doslim_name = doslimname_hash[doslim_id] if doslim_id in doslimname_hash else cancertype_hash[doslim_id]
                values = ["\""+doslim_name+"\""]
                for chr_id in chr_list + ["total"]:
                        val = count_hash["pass"][doslim_id][chr_id] if chr_id in count_hash["pass"][doslim_id] else 0
                        values.append(str(val))
                        count_hash["pass"][doslim_id]["total"] += val
                        if chr_id in count_hash["pass"]["total"]:
                                count_hash["pass"]["total"][chr_id] += val
                        else:
                                count_hash["pass"]["total"][chr_id] = val
                FW.write("%s\n" % (",".join(values)))
	FW.close()



if __name__ == '__main__':
        main()








