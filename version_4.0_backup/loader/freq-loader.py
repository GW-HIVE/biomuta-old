import os,sys
import string
import csv
import json
from optparse import Option_parser
import My_sQLdb



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
        config_json = json.loads(open(options.configfile, "r").read())
	progress_file = config_json["loader"]["outdir"] + "/freqloader-progress-"+src_type+".txt"


	DBH = My_sQLdb.connect(host = config_json["dbinfo"]["host"], 
				user = config_json["dbinfo"]["userid"], 
				passwd = config_json["dbinfo"]["password"],
				db = config_json["dbinfo"]["dbname"])
        cur = DBH.cursor()

	doslimid2cancerid = {}
	sql = "SELECT id,do_id FROM biomuta_cancer"
	cur.execute(sql)
	for row in cur.fetchall():
		doslimid2cancerid[row[1]] = row[0]

	chr_list = map(str,range(1,23)) + ["X","Y"]

	total_freq = 0
	n = 0
	report_progress(progress_file, "Started loading to biomuta_mutation_freq", True)	
	for chr_id in chr_list:
		comboid2mutid = {}
		sql = "SELECT id,chr,pos,ref,alt FROM biomuta_mutation where chr = '%s' " % (chr_id)
		cur.execute(sql)
		for row in cur.fetchall():	
			combo_id = ",".join([row[1],str(row[2]),row[3],row[4]])
			comboid2mutid[combo_id] = row[0]
		mutlist_file = config_json["datareporter"][src_type]["output_dir"]
		mutlist_file += "/mutlist.chr"+chr_id+".csv"
		if os.path.exists(mutlist_file) == False:
			continue
		with open(mutlist_file, 'rb') as tsvfile:
                      	tsvreader = csv.reader(tsvfile, delimiter=',', quotechar='"')
                      	row_count = 0
			for row in tsvreader:
				row_count += 1
				combo_id = ",".join(row[0:4])
				skip_list = [combo_id not in comboid2mutid, row[4] == "False"] 
				if True in skip_list:
					continue
				term_list = row[6].strip().split(" ")
				for term in term_list:
					doslim_id = term.split(":")[0]
					freq = int(term.split(":")[1])
					if int(doslim_id) >= 0:
						total_freq += freq
						n += 1
						cancer_id = doslimid2cancerid[doslim_id]
						mutation_id = comboid2mutid[combo_id]
						sql = "INSERT INTO biomuta_mutation_freq "
						sql += "(mutation_id, cancer_id, data_src, frequency) VALUES "
						sql += "(%s,%s,'%s',%s)"%(mutation_id,cancer_id,src_type,freq)
						cur.execute(sql)
						if n % 1000000 == 0:
							msg = "Loaded %s mutations" % (n)
							report_progress(progress_file, msg)
	report_progress(progress_file, "Total frequency count = %s" % (total_freq)) 
	DBH.commit()

if __name__ == '__main__':
        main()








