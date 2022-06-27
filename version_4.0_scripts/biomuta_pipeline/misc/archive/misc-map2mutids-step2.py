import os,sys
import string
import commands
from optparse import OptionParser
import glob
import json
import csv


__status__ = "Dev"
__version__ = "4.0"






##############################################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-i", "--configfile", action = "store", dest = "configfile", help = "Config file")
 

    (options,args) = parser.parse_args()
    for file in ([options.configfile]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    json_file = options.configfile
    config_json = json.loads(open(json_file).read())

    #chr_list = map(str,range(1,23)) + ["X","Y"]
    chr_list = ["1", "2"]
    for chr_id_out in chr_list:
        pattern = "outdir/mutid-*-chr%s.csv" % (chr_id_out)
        patient_list = {}
        count_dict = {}
        data_frame = {}
        src_list = []
        for in_file in glob.glob(pattern):
            file_name = in_file.split("/")[-1]
            src_type = file_name.split("-")[1]
            src_list.append(src_type)
            with open(in_file, 'r') as FR:
                row_count = 0
                for line in FR:
                    row_count += 1
                    row = line.strip().split(",")
                    combo_id = ",".join(row[:-3])
                    print combo_id

                    if combo_id not in data_frame:
                        data_frame[combo_id] = {}
                    data_frame[combo_id][src_type] = row[-1]


    sys.exit()


    src_list = sorted(set(src_list))
	row1 = ["chr_id","chr_pos","ref_nt","alt_nt","patient_id","uniprot_canonical_ac","aa_pos","ref_aa","alt_aa"]
	row1 += src_list

	row2 = ["uniprot_canonical_ac","aa_pos","ref_aa","alt_aa"]
	row2 += src_list

	out_file1 = "outdir/mutid-map.1.csv"
	out_file2 = "outdir/mutid-map.2.csv"

	FW1 = open(out_file1, "w")
	FW2 = open(out_file2, "w")

	FW1.write("%s\n" % (",".join(row1)))
	FW2.write("%s\n" % (",".join(row2)))

	for combo_id in data_frame:
		row1 = [combo_id]
		row2 = [",".join(combo_id.split(",")[4:])]
		for src_type in src_list:
			xref_id = data_frame[combo_id][src_type] if src_type in data_frame[combo_id] else ""
			row1.append(xref_id)
			row2.append(xref_id)
		FW1.write("%s\n" % (",".join(row1)))
		FW2.write("%s\n" % (",".join(row2)))

	FW1.close()
	FW2.close()


if __name__ == '__main__':
	main()


