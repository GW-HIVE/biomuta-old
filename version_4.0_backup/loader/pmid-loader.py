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

        (options,args) = parser.parse_args()
        for file in ([options.configfile]):
                if not (file):
                        parser.print_help()
                        sys.exit(0)

        config_json = json.loads(open(options.configfile, "r").read())
	pmidlist_file = config_json["loader"]["outdir"] + "/pmid.csv"


	DBH = My_sQLdb.connect(host = config_json["dbinfo"]["host"], 
				user = config_json["dbinfo"]["userid"], 
				passwd = config_json["dbinfo"]["password"],
				db = config_json["dbinfo"]["dbname"])
        cur = DBH.cursor()

	n = 0
	chr_list = map(str,range(1,23)) + ["X","Y"]
	for chr_id in chr_list:
		comboid2mutid = {}
		sql = "SELECT id,chr,pos,ref,alt FROM biomuta_mutation where chr = '%s' " % (chr_id)
		cur.execute(sql)
		for row in cur.fetchall():	
			combo_id = ",".join([row[1],str(row[2]),row[3],row[4]])
			comboid2mutid[combo_id] = row[0]
		with open(pmidlist_file, 'rb') as tsvfile:
                      	tsvreader = csv.reader(tsvfile, delimiter=',', quotechar='"')
                      	row_count = 0
			for row in tsvreader:
				row_count += 1
				combo_id = ",".join(row[0:4])
				skip_list = [combo_id not in comboid2mutid] 
				if True in skip_list:
					continue
				pm_id = row[4]
				mutation_id = comboid2mutid[combo_id]
				sql = "INSERT INTO biomuta_mutation_pmid "
				sql += "(mutation_id, pm_id) VALUES "
				sql += "(%s,'%s')"%(mutation_id,pm_id)
				#print sql
				cur.execute(sql)
				n += 1
		print "finished loading for ", chr_id
	DBH.commit()

if __name__ == '__main__':
        main()








