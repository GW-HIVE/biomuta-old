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



        src_type = "tcga"
        input_dir = config_json["datareporter"][src_type]["inputdir"]
        pattern = input_dir + "/*/*/*.vcf"
        file_list = glob.glob(pattern)
        tested_count = {}
        for f in file_list:
            cancer_type = f.split("/")[7]
            if cancer_type not in tested_count:
                tested_count[cancer_type] = 0
            tested_count[cancer_type] += 1
        

        chr_list = map(str,range(1,23)) + ["X","Y"]
        #chr_list = ["1"]

        for chr_id_out in chr_list:
            in_file = "outdir/mutid-tcga-chr%s.1.csv" % (chr_id_out)
            seen = {}
            positive_count = {}
            with open(in_file, 'r') as FR:
                reader = csv.reader(FR, delimiter=',', quotechar='"')
                row_count = 0
                for row in reader:
                    row_count += 1
                    combo_id = ",".join(row[0:5] + row[7:])
                    cancer_type = row[5]
                    patient_id = row[6]
                    
                    if combo_id not in seen:
                        seen[combo_id] = {}
                    if cancer_type not in seen[combo_id]:
                        seen[combo_id][cancer_type] = {}
                    if patient_id not in seen[combo_id][cancer_type]:
                        seen[combo_id][cancer_type][patient_id] = True
                        if combo_id not in positive_count:
                            positive_count[combo_id] = {}
                        if cancer_type not in positive_count[combo_id]:
                            positive_count[combo_id][cancer_type] = 0
                        else:
                            positive_count[combo_id][cancer_type] += 1

            out_file = "outdir/mutid-tcga-chr%s.2.csv" % (chr_id_out)
            FW = open(out_file, "w")
            for combo_id in positive_count:
                for cancer_type in positive_count[combo_id]:
                    patients_positive = positive_count[combo_id][cancer_type]
                    patients_tested = tested_count[cancer_type]
                    FW.write("%s,%s,%s,%s\n" % (combo_id,cancer_type,patients_positive,patients_tested))
            FW.close()




if __name__ == '__main__':
	main()


