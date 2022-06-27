import os,sys
import string
import csv
import json
from optparse import Option_parser
import glob



__version__="1.0"
__status__ = "Dev"




###############################
def main():

	file_list = glob.glob("*.py")
	for file in file_list:
		with open(file, "r") as FR:
			for line in FR:
				line = line.replace("\"", "")
				line = line.replace("\,", "")
				line = line.replace(";", "; ")
				line = line.replace("[", "[ ")
				line = line.replace("]", "] ")
				line = line.replace("(", "( ")
				line = line.replace(")", ") ")
				line = line.replace(".", ". ")
				word_list = line.strip().split(" ")
				for word in word_list:
					for i in xrange(0, len(word)-1):
						if word[i].lower() == word[i] and word[i+1].lower() != word[i+1]:
							print word.strip()





if __name__ == '__main__':
        main()








