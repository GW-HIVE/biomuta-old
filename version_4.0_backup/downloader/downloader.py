import json
import commands
import sys
import paramiko
import os
import csv
from glob import glob

__version__="1.0"
__status__ = "Dev"


###########################

#Data Resource
class DataResource(object):

        def __init__(self, config_file_path):
                self.config_json = json.loads(open(config_file_path).read())["downloader"]

        def download(self, source):

                # SOURCE TYPE FTP or HTTPS
                if (self.config_json[source]["sourceprotocol"] == "ftp" or self.config_json[source]["sourceprotocol"] == "https"):
                        progress_file = self.config_json[source]["outdir"] + "progress.txt"
                        for i in range(len(self.config_json[source]["sourceurilist"])):
                                download_url = self.config_json[source]["sourceurilist"][i]
                                output_file = self.config_json[source]["downloadfolder"] + "/" + os.path.basename(download_url)
                                log_file = self.config_json[source]["outdir"] + "/logfile." + str(i) + ".txt"
                                FW = open(progress_file, "a")
                                FW.write("Downloading file %s\n" % (i))
                                cmd = "wget -O " + output_file + " " + download_url + " >& " + log_file
                                x = commands.getoutput(cmd)
                                FW.write("Downloaded File %s\n" % (i))
                                FW.write("Unzipping file %s\n" % (i))
                                cmd = "gunzip " + output_file
                                x = commands.getoutput(cmd)
                                FW.write("Unzipped File %s\n" % (i))
                                FW.close()
                        print "Download Complete"

                # SOURCE TYPE SFTP
                elif self.config_json[source]["sourceprotocol"] == "sftp":
                        host = self.config_json[source]["host"]
                        user_login = self.config_json[source]["userlogin"]
                        user_pass = self.config_json[source]["userpass"]
                        port = 22
                        transport = paramiko.Transport((host, port))
                        transport.connect(username=user_ogin, password=user_pass)
                        sftp = paramiko.SFTPClient.from_transport(transport)
                        progress_file = self.config_json[source]["outdir"] + "progress.txt"

                        for i in range(len(self.config_json[source]["sourceurilist"])):
                                FW = open(progress_file, "a")
                                FW.write("Downloading file %s\n" % (i))
                                file_path = self.config_json[source]["sourceurilist"][i]
                                local_path = self.config_json[source]["downloadfolder"] + os.path.basename(file_path)
                                sftp.get(file_path, local_path)
                                FW.write("Downloaded File %s\n" % (i))
                                FW.write("Unzipping file %s\n" % (i))
                                cmd = "gunzip " + local_path
                                x = commands.getoutput(cmd)
                                FW.write("Unzipped File %s\n" % (i))
                                FW.close()

                        sftp.close()
                        transport.close()
                        print "Download Complete"

                # SOURCE TYPE MANIFEST
                elif self.config_json[source]["sourceprotocol"] == "manifest":

                        project_list_file = self.config_json[source]["project_list_file"]
                        token_file = self.config_json[source]["token_file"]
                        download_prg = self.config_json[source]["download_prg"]
                        download_dir = self.config_json[source]["downloadfolder"]
                        FR = open(project_list_file, "r")

                        for line in FR:
                                prj_name = line.strip()
                                if prj_name != "" and prj_name[0:1] != "#":
                                        manifest_file = download_dir + "manifest-files/" + prj_name + ".reduced.tsv"
                                        output_dir = download_dir + prj_name + "/"
                                        cmd = download_prg + " download -t " + token_file + " -m " + manifest_file
                                        cmd += " -d " + output_dir
                                        x = commands.getoutput(cmd)
                        FR.close()
                        print "Download Complete"

                # SOURCE TYPE MANUAL
                elif self.config_json[source]["sourceprotocol"] == "manual":

                        print "\n\n"
                        print "\tThis source has not yet been automated\n"
                        print "\tTo Manually Download:\n\n"
                        print "\t\t----------------------------------------------------"
                        print "\t\tDownload URL : ", self.config_json[source]["sourceurilist"][0]
                        print "\t\t----------------------------------------------------"
                        print "\t\tUser Name : ", self.config_json[source]["userlogin"]
                        print "\t\t----------------------------------------------------"
                        print "\t\tUser Pass : ", self.config_json[source]["userpass"]
                        print "\t\t----------------------------------------------------\n\n"
                        print "\tSave the files in the following location:\n"
                        print "\t", self.config_json[source]["downloadfolder"]
                        print "\t\n\n"

        def checkmd5(self, source):
                digest_dict = {}
                file_list = []
                checksum_results = self.config_json[source]["outdir"] + "checksum.txt"

                manifest_pattern, vcf_pattern = "", ""
                if source == "tcga":
                        manifest_pattern = self.config_json[source]["downloadfolder"] + "/manifest-files/*.reduced.tsv"
                        vcf_pattern = self.config_json[source]["downloadfolder"] + "*/*/*.vcf"

                digest_dict = {}
                for f in glob(manifest_pattern):
                        with open(f, 'rb') as tsvfile:
                                tsvreader = csv.reader(tsvfile, delimiter='\t',quotechar='"')
                                row_count = 0
                                for row in tsvreader:
                                        row_count +=1
                                        if row_count ==1:
                                                continue
                                        file_name, md5sumref = row[1], row[2]
                                        digest_dict[file_name] = md5sumref

		header_list = ["flag","filename","calculated_md5","reference_md5"]
		FW = open(checksum_results, "w")
		FW.write("%s\n" % (",".join(header_list)))
                for f in glob(vcf_pattern):
                        cmd = self.config_json[source]["md5Prg"] + " " + f
                        checkmd5sum = commands.getoutput(cmd).split(" ")[0].strip()
                        check_file_name = os.path.basename(f)
                	flag = checkmd5sum == digest_dict[check_file_name]
			FW.write("%s,%s,%s,%s\n" % (flag, f, checkmd5sum, digest_dict[check_file_name]))
                FW.close()
		


