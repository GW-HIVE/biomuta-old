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

	chr_list = map(str,range(1,23)) + ["X","Y"]
        #chr_list = ["1"]


	coding_dict = {}
	print "Started loading coding_dict dict ... "
	mutlist_file = config_json["mutmapper"]["outdir"] + "/mutlist-merged.csv"
	with open(mutlist_file, 'r') as FR:
                reader = csv.reader(FR, delimiter=',', quotechar='"')
                row_count = 0
                for row in reader:
                        if row[-1] == "coding":
				mut_index,chrom,pos,ref,alt = row[0:5]
				if chrom not in chr_list:
                                    continue
                                combo_id = "%s,%s,%s,%s" % (chrom,pos,ref,alt)
				coding_dict[combo_id] = True
	print "done!"


	eff_dict = {}
	print "Started loading eff_dict ... "
	seen = {"aacombo":{}}
	mutmap_file =  config_json["mutmapper"]["outdir"] + "/mutmap.2.csv"
        with open(mutmap_file, 'r') as FR:
                reader = csv.reader(FR, delimiter=',', quotechar='"')
                row_count = 0
                for row in reader:
			row_count += 1
			mut_index,chrom,pos,ref,alt = row[0:5]
                        if chrom not in chr_list:
                            continue
			refaa,altaa = row[12], row[13]
			ac,pos_aa, is_canon = row[14], row[15], row[20]
			combo_id = "%s,%s,%s,%s" % (chrom,pos,ref,alt)
			if is_canon == "True" and combo_id in coding_dict:
				if combo_id not in eff_dict:
					eff_dict[combo_id] = []	
				aa_combo_id = "%s,%s,%s,%s,%s" % (combo_id,ac,pos_aa,refaa,altaa)
				if aa_combo_id not in seen["aacombo"]:
					eff_dict[combo_id].append({"ac":ac, "pos":pos_aa, "refaa":refaa, "altaa":altaa})
					seen["aacombo"][aa_combo_id] = True
	print "done!"





        src_type = "icgc"
	input_dir = config_json["datareporter"][src_type]["inputdir"]
        pattern = input_dir + "/*.vcf"
        file_list = glob.glob(pattern)
	
        for chr_id_out in chr_list:
		file_count = 0
		mutid_dict = {}
		print "started %s of %s ... " % (chr_id_out, src_type)

                file_name = "mutid-%s-chr%s.1.csv" % (src_type, chr_id_out)
                out_file =  "outdir/" + file_name
                FW = open(out_file, "w")
		for vcf_file in file_list:
			file_count += 1
			with open(vcf_file, 'r') as FR:
               			row_count = 0
				for line in FR:
					row_count += 1
					row = line.strip().split("\t")
					if row[0][0:1] == "#":	
						continue
					chr_id = row[0]
					if chr_id != chr_id_out:
						continue
                                        mut_id = ""
					if row[2] not in ["", "."]:
                                            mut_id = row[2]
                                        att_dict = {}
                                        for att in row[-1].split(";"):
                                            parts = att.split("=")
                                            if len(parts) == 2:
                                                att_dict[parts[0]] = parts[1]
                                        
                                        nt_row = [row[0],row[1],row[3],row[4],mut_id]
                                        combo_id = ",".join(nt_row[0:4])
                                        if combo_id in eff_dict:
                                            for obj in eff_dict[combo_id]:
                                                aa_row = [obj["ac"], obj["pos"],obj["refaa"],obj["altaa"]]
                                                for freq_info in  att_dict["OCCURRENCE"].split(","):
                                                    cancer_type = freq_info.split("-")[0]
                                                    n_val = freq_info.split("|")[1]
                                                    d_val = freq_info.split("|")[2]
                                                    new_row = nt_row + aa_row + [cancer_type,n_val,d_val]
                                                    FW.write("%s\n" % (",".join(new_row)))
		FW.close()
		print "done!"



if __name__ == '__main__':
	main()


