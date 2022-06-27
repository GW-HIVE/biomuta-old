import os,sys
import string
import csv
import json
from optparse import OptionParser
from Bio import SeqIO
from Bio.Seq import Seq


__version__="1.0"
__status__ = "Dev"








###############################
def main():

	usage = "\n%prog  [options]"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--configfile",action="store",dest="configfile",help="NT file")


	(options,args) = parser.parse_args()
	for file in ([options.configfile]):
		if not (file):
			parser.print_help()
			sys.exit(0)

	config_json = json.loads(open(options.configfile, "r").read())
        grch_ver = config_json["grchversion"]
        uniprot_ver = config_json["uniprotversion"]

        idmapping_file = config_json["ref"][uniprot_ver]["idmapping"]
        

        pepid2isoformac_file = "outdir/pepid2isoformac.csv"
        ranges_file = "outdir/map2uniprot.ranges.csv"
        dbsnp_in_file_one = "outdir/dbsnp.1.csv"
        dbsnp_in_file_two = "outdir/dbsnp.2.csv"
        dbsnp_out_file = "outdir/dbsnp.3.csv"



	is_canon = {}
	with open(idmapping_file, "r") as FR:
            data_grid = csv.reader(FR, delimiter=',', quotechar='"')
            row_count = 0
            for row in data_grid:
                isoform = row[0]
                is_canon[isoform] = True 




	pepid2acmap = {"uniprot1":{}, "uniprot2":{}}
	with open(pepid2isoformac_file, "r") as FR:
            data_grid = csv.reader(FR, delimiter=',', quotechar='"')
            row_count = 0
            for row in data_grid:
                pep_id, ac = row[0], row[2]
                if row[1] == "exact":
                    if pep_id not in pepid2acmap["uniprot1"]:
                        pepid2acmap["uniprot1"][pep_id] = []
                    pepid2acmap["uniprot1"][pep_id].append(ac)



	with open(ranges_file, "r") as FR:
                data_frame = csv.reader(FR, delimiter=',', quotechar='"')
                row_count = 0
		for row in data_frame:
                        row_count += 1
			if row_count == 1:
				continue
			pep_id, ac = row[0], row[1]
			pepid2acmap["uniprot2"][pep_id] = {}
			for r in row[2].strip().split("|"):
				r1, r2 = r.split(":")
				start_pos1,end_pos1 = r1.split("-")
				start_pos2,end_pos2 = r2.split("-")
				n = 0
				for pos in xrange(int(start_pos1), int(end_pos1) + 1):
					val = "%s %s" % (ac,int(start_pos2) + n)
					pepid2acmap["uniprot2"][pep_id][str(pos)] = val
					n += 1


	map_dict = {}
	with open(dbsnp_in_file_two, 'rb') as FR:
		data_frame = csv.reader(FR, delimiter=',', quotechar='"')
		row_count = 0
		for row in data_frame:
                        row_count += 1
			if row_count == 1:
				continue
			pep_id = row[6]
			position_in_pep = row[8]
			uniprot_ac, pos_in_uniprot, flag = "-", "-","-"
			if pep_id in pepid2acmap["uniprot2"]:
				if position_in_pep in pepid2acmap["uniprot2"][pep_id]:
					uniprot_ac = pepid2acmap["uniprot2"][pep_id][position_in_pep].split(" ")[0]
					pos_in_uniprot = pepid2acmap["uniprot2"][pep_id][position_in_pep].split(" ")[1]
					flag = "nonexact"
			if pep_id in pepid2acmap["uniprot1"]:
				uniprot_ac = pepid2acmap["uniprot1"][pep_id][0]
				pos_in_uniprot = position_in_pep
				flag = "exact"
			newrow = row[0:6] + row[7:8] + row[6:7] + row[8:9] + [uniprot_ac, pos_in_uniprot,row[-2],row[-1]]
			combo_id = ",".join(newrow[1:5])
			map_dict[combo_id] = newrow[5:]



	field_list = ["chr", "pos", "id", "ref", "alt", "caf", "sao", "geneinfo","vc", "nsf","nsm", "nsn", "syn"]
	field_list += ["trsid","pos_in_trs", "pep_id", "pos_in_pep", "uniprot_isoform_ac", "pos_in_uniprot_isoform", "ref_aa", "alt_aa"]
	FW = open(dbsnp_out_file, "w")
	FW.write("\"%s\"\n" % ("\",\"".join(field_list)))
	with open(dbsnp_in_file_one, 'rb') as tsvfile:
                tsvreader = csv.reader(tsvfile, delimiter=',', quotechar='"')
                row_count = 0
                for row in tsvreader:
                        row_count += 1
                        if row_count == 1:
                                continue
			combo_id = "%s,%s,%s,%s" % (row[0], row[1],row[3], row[4])
			if combo_id in map_dict:
				newrow = row + map_dict[combo_id]
				FW.write("\"%s\"\n" % ("\",\"".join(newrow)))
			if row_count % 1000000 == 0:
                                print row_count

	FW.close()

if __name__ == '__main__':
        main()








