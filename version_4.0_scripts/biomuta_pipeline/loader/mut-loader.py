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
	parser.add_option("-i","--configfile",action="store",dest="configfile",help="NT file")


	(options,args) = parser.parse_args()
	for file in ([options.configfile]):
		if not (file):
			parser.print_help()
			sys.exit(0)




	config_json = json.loads(open(options.configfile, "r").read())
 	DBH = My_sQLdb.connect(host = config_json["dbinfo"]["host"], 
				user = config_json["dbinfo"]["userid"], 
				passwd = config_json["dbinfo"]["password"],
				db = config_json["dbinfo"]["dbname"])
        cur = DBH.cursor()


	seen = {}

	mutlist_file = config_json["mutmapper"]["outdir"] + "/mutlist-merged.csv"
	progress_file = config_json["loader"]["outdir"] + "/mutloader-progress.txt"

	report_progress(progress_file, "Started loading to biomuta_mutation", True)
	#Insert placeholder first
	sql = "INSERT INTO biomuta_mutation (chr, pos, ref, alt) VALUES ('-',0,'-','-') "
	cur.execute(sql)
	with open(mutlist_file, 'rb') as tsvfile:
                tsvreader = csv.reader(tsvfile, delimiter=',', quotechar='"')
                row_count, n = 0, 0
                for row in tsvreader:
                        row_count += 1
                        if row_count == 1:
                                continue
			if row[-1] == "noncoding":
				continue

			mut_id, chr_id, pos, ref, alt = row[0], row[1],row[2], row[3], row[4]
			key = "%s,%s,%s,%s" % (chr_id.strip(), pos.strip(), ref.strip(), alt.strip())
			if key in seen:
				continue
			seen[key] = True
			n += 1
			sql = ("INSERT INTO biomuta_mutation (chr, pos, ref, alt)  VALUES ('%s',%s,'%s','%s') " % (chr_id, int(pos), ref, alt))
			cur.execute(sql)
			if n % 1000000 == 0:
				report_progress(progress_file, "Loaded %s mutations" % (n))

	DBH.commit()

if __name__ == '__main__':
        main()








