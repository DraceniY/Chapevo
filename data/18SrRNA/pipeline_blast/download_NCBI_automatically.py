# This script parse eukaryotes.csv and match names in this file with list of eukaryote
# download all transcript
import sys
import pandas as pd
import numpy as np
import glob , os , re
from ftplib import FTP
import subprocess
from subprocess import call

ncbi_file = "/home/yayu/Documents/Maitrise_Bio-informatique/paper/data/eukaryotes.csv"
liste_eukaryote = "/home/yayu/Documents/Maitrise_Bio-informatique/paper/data/HSP_FINAL.csv"

def parse_files ( liste_eukaryote , ncbi_file ) :
	liste_euka=[]
	eukaryote_NCBI=[]
	euka_download=[]
	with open ( liste_eukaryote , "r") as f :
		line = f.readline()
		for line in f :
			ID_euka = line.split(",") [1]
			liste_euka.append(ID_euka)

	with open ( ncbi_file , "r") as f1:
		line1 = f1.readline()
		for line1 in f1:
			ID_euka_NCBI1 = ((line1.split(",") [0]).split(" ")[0]).replace('"' , '')
			ID_euka_NCBI2 = ((line1.split(",") [0]).split(" ")[1]).replace('"' , "")
			adress_web =  line1.split (",")[15]
			eukaryote_NCBI.append("_".join([ID_euka_NCBI1,ID_euka_NCBI2]))
			euka_download.append(adress_web.replace("\n" , "").replace("''" , ""))

# Search match ID eukaryote
	ID_not_found = []
	ID_find = {}

	for i in liste_euka:
		if (i in eukaryote_NCBI) and (i not in ID_find.keys()) :
			if (str(euka_download[eukaryote_NCBI.index(i)]) == '\n'):
				pass
			else:
				ID_find[i] = [euka_download[eukaryote_NCBI.index(i)]]
		if i not in eukaryote_NCBI :
			ID_not_found.append(i)

	not_have_transcript = []
	# ftp = FTP('ftp.ncbi.nih.gov')   # Access to NCBI ftp through Python script
	# ftp.login()

	for i , j in ID_find.items():
		try :
			url = j[0].split(".gov/")[1].replace('"' , '')
			fetch = "GCF_"+j[0].split("/GCF_")[1].replace('"' , '')
			subprocess.call('curl -O https://ftp.ncbi.nih.gov/'+url+"/"+fetch+'_rna.fna.gz' , shell=True)
			# subprocess.check_output('curl -O https://ftp.ncbi.nih.gov/'+url+"/"+fetch+'_rna.fna.gz')

		except :
			print i
			not_have_transcript.append(i)
			# ftp.cwd(url)
			# for currentfile in ftp.nlst():
			# 	if fetch+"_rna.fna.gz" in currentfile:
			# 		print("Found : "+i)
			# 		print('curl -O https://ftp.ncbi.nih.gov/'+url+"/"+fetch+'_rna.fna.gz')
			# 	else:
			# 		not_have_transcript.append(i)

	os.system("gunzip *.gz")
	os.system("mkdir download")
	os.system ("mv *.fna download")

parse_files(liste_eukaryote , ncbi_file)
