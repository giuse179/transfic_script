#!/soft/bin/python
# COMANDI DA SHELL
#sort -k1,1 -k2,2n -k3,3n snp10.vep.bed | bgzip > snp10.vep.bed.gz
#home/giumar/soft/tabix-0.2.6/tabix -p bed snp10.vep.bed.gz
#perl /home/giumar/soft/variant_effect_predictor/variant_effect_predictor.pl -i prova.bed --cache --offline --chr 1 --protein --sift=s --polyphen=s -o prova_VEP
import sys, csv, re
from os import listdir, system
from os.path import isfile, join

path_soft_VEP='variant_effect_predictor.pl' ### add choose to path
path_VEP='/variant_effect_predictor/.vep'
path_soft_transfic='transf_scores.pl' ### add choose to path
path_soft_fathmm='fathmm.py' ### add choose to path

def run_VEP_1000genome(path_soft,path_VEP):
  path_dataset='/transfic_project/1000genomes/1000_bedfile_VEP/'  ### add argv choose
	path_outfile='/transfic_project/1000genomes/VEP_result/'
	list_file=[f for f in listdir(path_dataset) if isfile(join(path_dataset,f))]	
	for data in list_file:
		z=re.compile(r"\d+")
		a=re.search(z,data)
		if a == None:k='X'
		else: k=a.group() 						
		if k =='X' or k=='12' or k=='13': 
			print data,'......start!!!'					
			system ("perl %s -i %s%s --cache --offline --chr %s --protein --sift=b --polyphen=b --buffer_size 20000 --dir %s -o %s%s_VEP" % (path_soft,path_dataset,data,k,path_VEP,path_outfile,data))
			print data,'....finished!!!'
			print '--------------------------------------------------------------------------------'
		else: continue

def transfic_process_cosmic(path_soft):
	path_dataset='/transfic_project/file_per_transfic/cosmic/' ### add argv choose
	path_outfile='/transfic_project/transfic_result/cosmic/'
	list_file=[f for f in listdir(path_dataset) if isfile(join(path_dataset,f))]
	for data in list_file:
		print data,'......start!!!'		
		system ("perl %s gosmf %s%s %s%s_tranfic" % (path_soft,path_dataset,data,path_outfile,data))
		print data,'   finished!!!'

def fathmm_process_cosmic(path_soft):
	path_dataset='/transfic_project/file_per_fathHMM/cosmic/' ### add argv choose
	path_outfile='/transfic_project/fathHMM_result/cosmic/'
	list_file=[f for f in listdir(path_dataset) if isfile(join(path_dataset,f))]
	for data in list_file:
		print data,'......start!!!'		
		system ("/soft/bin/python  %s -i %s%s -o %s%s" % (path_soft,path_dataset,data,path_outfile,data))
		print data,'....finished!!!'

def fathmm_process_VEP(path_soft):
	path_dataset='/transfic_project/file_per_fathHMM/VEP/'
	path_outfile='/transfic_project/fathHMM_result/VEP/'
	list_file=[f for f in listdir(path_dataset) if isfile(join(path_dataset,f))]
	for data in list_file:
		print data,'......start!!!'		
		system ("/soft/bin/python  %s -i %s%s -o %s%s" % (path_soft,path_dataset,data,path_outfile,data))
		print data,'....finished!!!'

if __name__=="__main__":
#	run_VEP_1000genome(path_soft_VEP,path_VEP)
#	transfic_process_cosmic(path_soft_transfic)
#	fathmm_process_cosmic(path_soft_fathmm)
#	fathmm_process_VEP(path_soft_fathmm)

