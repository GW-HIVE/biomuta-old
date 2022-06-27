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
 
	map_file = config_json["mutmapper"]["outdir"] + "/mutmap.2.csv"
	progress_file = config_json["loader"]["outdir"] + "/effloader-progress.txt"


	DBH = My_sQLdb.connect(host = config_json["dbinfo"]["host"], 
				user = config_json["dbinfo"]["userid"], 
				passwd = config_json["dbinfo"]["password"],
				db = config_json["dbinfo"]["dbname"])
        cur = DBH.cursor()


	chr_list = map(str,range(1,23)) + ["X","Y"]

	report_progress(progress_file, "Started loading to biomuta_mutation_eff", True)	
	cumm = 0
	for chr_id in chr_list:
		comboid2mutid = {}
		n = 0
		sql = "SELECT id,chr,pos,ref,alt FROM biomuta_mutation where chr = '%s' " % (chr_id)
		cur.execute(sql)
		for row in cur.fetchall():	
			combo_id = ",".join([row[1],str(row[2]),row[3],row[4]])
			comboid2mutid[combo_id] = row[0]
		with open(map_file, 'rb') as tsvfile:
                     	tsvreader = csv.reader(tsvfile, delimiter=',', quotechar='"')
                     	row_count = 0
			for row in tsvreader:
                              	row_count += 1
				combo_id = ",".join(row[1:5])
				skip_list = [row[1] != chr_id, combo_id not in comboid2mutid] 
				if True in skip_list:
					continue
				n += 1
				mutation_id = comboid2mutid[combo_id]
				transcsript_id = row[5]
				peptide_id = row[6]
				pos_in_cds = row[7] 
				pos_in_pep = row[8]
				pos_in_codon = row[9]
				ref_codon = row[10]
				alt_codon = row[11]
				ref_residue = row[12]
				alt_residue = row[13]
				canonical_ac = row[14]
				pos_in_canonical_ac = row[15] if row[15] != "-" else -1
				flag1 = row[16]
				refseq_ac = row[17]
				pos_in_ref_ac = row[18] if row[18] != "-" else -1
				flag2 = row[19]
				is_canonical = row[20]
				if is_canonical == "False":
					continue
				
				sql = "INSERT INTO biomuta_mutation_eff "
				sql += "(mutation_id,transcript_id,peptide_id,pos_in_cds,pos_in_pep,pos_in_codon,ref_codon,alt_codon,ref_residue,alt_residue,canonical_ac,pos_in_canonical_ac,refseq_ac,pos_in_ref_ac) VALUES "
				sql += "(%s,'%s','%s',%s,%s,%s,'%s','%s','%s','%s','%s',%s,'%s',%s)" % (mutation_id,transcsript_id,peptide_id,pos_in_cds,pos_in_pep,pos_in_codon,ref_codon,alt_codon,ref_residue,alt_residue,canonical_ac,pos_in_canonical_ac,refseq_ac,pos_in_ref_ac)
				cur.execute(sql)
				#print sql + "\n"

		cumm += n
		msg = "%s: %s mutations (%s so far)" % (chr_id, n, cumm)
		report_progress(progress_file, msg)
	DBH.commit()

if __name__ == '__main__':
        main()








