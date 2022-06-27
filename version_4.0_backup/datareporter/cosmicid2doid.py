import csv
import json
import sys
from optparse import OptionParser

__version__ = "4.0"
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



############################################################3333
def main():
	
	usage = "\n%prog [options]"
	parser = OptionParser(usage, version = "%prog" + __version__)
	parser.add_option("-i", "--configfile", action = "store", dest = "configfile" , help = "Config Json File Path")
	
	(options, args) = parser.parse_args()
	for file in ([options.configfile]):
		if not(file):
			parser.print_help()
			sys.exit(0) 

	config_file = options.configfile
	config_json = json.loads(open(config_file).read())
	
	in_file_cosmic1 = config_json["datareporter"]["cosmic"]["inputdir"] + "CosmicMutantExport.tsv"
	in_file_cosmic2 = config_json["datareporter"]["cosmic"]["inputdir"] + "CosmicNCV.tsv"
	out_file1 = config_json["datareporter"]["cosmic"]["outputdir"] + "cosmic-mutid2doid.tsv"
	out_file2 = config_json["datareporter"]["cosmic"]["outputdir"] + "cosmic-withoutdoidslim.tsv"
	progress_file = config_json["datareporter"]["cosmic"]["outputdir"] + "mutid2doid-progress.tsv"
	doid_file = config_json["datareporter"]["doid"]["outputdir"] + "/doid-mapping.reduced.tsv"
        
	
	seen = {"cosmid":{}, "cancertype1":{}, "doid":{}, "cancertype2":{}}

	#Parsing from the Comic tsv - Mapping Cosmic_id to Subtype for Coding 
	cosmid2cancertypelist = {}
	report_progress(progress_file, "Mapping Cosmic_id to Subtype for Coding", True)
	with open(in_file_cosmic1, 'r') as FR:
                csv_grid = csv.reader(FR, delimiter='\t', quotechar='"')
                row_count = 0
		for row in csv_grid:
			row_count +=1
			if row_count == 1 :
				continue
			cosm_id = row[16]
			seen["cosmid"][cosm_id] = True
			cancer_type_list = row[7:15]
			cosmid2cancertypelist[cosm_id] = cancer_type_list
			for cancer_type in cancer_type_list:
				seen["cancertype1"][cancer_type] = True
			if row_count % 1000000 == 0:
				report_progress(progress_file, "Processed %s rows" % (row_count))

	#Parsing from the Comic tsv - Mapping Cosmic_id to Subtype for Non_coding
	with open(in_file_cosmic2, 'r') as FR:
                csv_grid = csv.reader(FR, delimiter='\t', quotechar='"')
                row_count = 0
                for row in csv_grid:
                        row_count +=1
                        if row_count == 1 :
                                continue
                        cosm_id = row[11]
                        seen["cosmid"][cosm_id] = True
                        cancer_type_list = row[3:11]
                        cosmid2cancertypelist[cosm_id] = cancer_type_list 
                        for cancer_type in cancer_type_list:
                                seen["cancertype1"][cancer_type] = True
                        if row_count % 1000000 == 0:
                                report_progress(progress_file, "Processed %s rows" % (row_count))

	

	 
	#Load DOID mappings
	doslimid_hash = {}
	doslimname_hash = {}
        load_do_slim (doid_file, doslimid_hash, doslimname_hash, [])
        #Repeat loading to overwrite mapping to 305 if more specific mapping exists
        load_do_slim (doid_file, doslimid_hash, doslimname_hash, ["305"])


	report_progress(progress_file, "Unique DOIDs seen = ",  len(seen["doid"].keys()))
        report_progress(progress_file, "Unique subtypes seen in second file = ",  len(seen["cancertype2"].keys()))


	FW1 = open(out_file1, "w")
	FW2 = open(out_file2, "w")
	#Mapping Do_id to Comic_id Using Subtype
	status = {}
	report_progress(progress_file, "Mapping Do_id to Comic_id Using Subtype")
	for cosm_id in cosmid2cancertypelist:
               	flag = False
		for cancer_type in cosmid2cancertypelist[cosm_id]:
                       	if cancer_type in doslimid_hash:
				do_id = doslimid_hash[cancer_type]
				flag = True
				FW1.write("%s\t%s\t%s\n" % (cosm_id,do_id,cancer_type))
				break
		if flag == False:
			FW2.write("%s\t%s\t%s\n" % (cosm_id,-1,";".join(cosmid2cancertypelist[cosm_id])))
	FW1.close()
	FW2.close()			
	
	

if __name__ == "__main__":
	main()



