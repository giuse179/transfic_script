#!/soft/bin/python
#12/04/2013

import sys, csv, re
import pandas as pd
from os import listdir, system
from os.path import isfile, join
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

path_soft_VEP='/genomics/users/giumar/soft/variant_effect_predictor/variant_effect_predictor.pl' ### add choose to path
path_VEP='/genomics/users/giumar/soft/variant_effect_predictor/.vep'
path_soft_transfic='/home/giumar/soft/transfic/bin/transf_scores.pl' ### add choose to path
path_soft_fathmm='/home/giumar/soft/fathmm-master/cgi-bin/fathmm.py' ### add choose to path

#########################################################################################################################
def bedfile_transform():
  path_project='/genomics/users/giumar/transfic_project/1000genomes/'
	list_file=[f for f in listdir(path_project+'original_data/') if isfile(join(path_project+'original_data/',f))]		
	for dataset in list_file:		
		filein=path_project+'1000genome_originalDATA/'+dataset		
		print filein,'...start!!!'		
		filename_out= path_project+'1000genome_bedfile/'+'%s_mut' %(dataset)
		fileout=open(filename_out,'w')		
		dataframes=pd.read_table(filein,sep='\t',names=['chr','pos','wt','mt'])
		dataframes['mut']=['/'.join(row) for row in dataframes[dataframes.columns[2:]].values]		
		dataframes['strand']= '+'	
		#print dataframes[:5]		
		cols=['chr','pos','pos','mut','strand']
		changed_dataframe=dataframes.ix[:,cols]		
		changed_dataframe.to_csv(fileout,sep='\t',index=False,header=False)		
		fileout.close()			
		print fileout,('...finish!!!')
########################################################################################################################
def run_VEP_1000genome(path_soft,database_VEP): ### bottlenek
	path_dataset='/genomics/users/giumar/transfic_project/1000genomes/1000_bedfile_VEP/'  ### add argv choose
	path_outfile='/genomics/users/giumar/transfic_project/1000genomes/VEP_result/'
	list_file=[f for f in listdir(path_dataset) if isfile(join(path_dataset,f))]	
	for data in list_file:
		z=re.compile(r"\d+")
		a=re.search(z,data)
		if a == None:k='X'
		else: k=a.group() 						
		if k =='X' or k=='12' or k=='13': #### aggiunta condizione per dividere il processo per CROMOSOMA in modo da lanciarlo in piu PC
			print data,'......start!!!'					
			system ("perl %s -i %s%s --cache --offline --chr %s --protein --sift=b --polyphen=b --buffer_size 20000 --dir %s -o %s%s_VEP" % (path_soft,path_dataset,data,k,database_VEP,path_outfile,data))
			print data,'....finished!!!'
			print '--------------------------------------------------------------------------------'
		else: continue
###########################################################################################################################
def from_VEP_to_FATHMM_chunk():	
	path_dataset='/genomics/users/giumar/transfic_project/1000genomes/VEP_result/'  ### add argv choose??????
	path_outfile='/genomics/users/giumar/transfic_project/file_per_fathHMM/VEP/'
	list_file=[f for f in listdir(path_dataset) if isfile(join(path_dataset,f))]	
	#list_file=['snp1_mut_VEP']
	for dataset in list_file:	
		filein=path_dataset+dataset
		filename_out=path_outfile+'%s_fathmm2' % (dataset)	
		print dataset,'...start!!!'
		fileout=open(filename_out,'w')	
		cols=['Uploaded_variation','Location','Allele','Gene','Feature','Feature_type','Consequence','cDNA_position','CDS_position','Protein_position','Amino_acids','Codons','Existing_variation','Extra']		
		chunk_dataframes=pd.read_csv(filein,sep='\t',header=None,na_values=['-'],index_col=[6],skiprows=12,names=cols,chunksize=100000)
		i=0		
		for dataframes in chunk_dataframes:			
			cols2=['Protein_position','Amino_acids','Extra']
			dataframes=dataframes.ix[:,cols2].copy()	
			select=dataframes.index=='missense_variant'			
			dataframes=dataframes.ix[select]				
			dataframes['wild_Amic']=dataframes['Amino_acids'].str.split('/').str.get(0)
			dataframes['mut_Amic']=dataframes['Amino_acids'].str.split('/').str.get(1)			
			a=dataframes['Extra'].str.split(';').str.get(0)	
			dataframes['prot_name']=a.str.split('=').str.get(1)	
			dataframes['intpos']=np.array(dataframes.Protein_position).astype(np.int32)	
			dataframes['pos']=np.array(dataframes.intpos).astype('|S5')
			cols3=['prot_name','wild_Amic','pos','mut_Amic']			
			dataframes=dataframes.ix[:,cols3].copy()
			dataframes['mut']=[''.join(filter(None,row)) for row in dataframes[dataframes.columns[1:]].values]	
			dataframes.to_csv(fileout,sep='\t',index=False,header=False,cols=['prot_name','mut'])		
			i+=1			
			print dataset,'...chunk--->',i
		print fileout,('.........finish!!!')
		fileout.close()
##############################################################################################################################
def fathmm_process_VEP(path_soft):
	path_dataset='/genomics/users/giumar/transfic_project/file_per_fathHMM/VEP/'
	path_outfile='/genomics/users/giumar/transfic_project/fathHMM_result/VEP/'
	list_file=[f for f in listdir(path_dataset) if isfile(join(path_dataset,f))]
	for data in list_file:
		print data,'......start!!!'		
		system ("/soft/bin/python  %s -i %s%s -o %s%s" % (path_soft,path_dataset,data,path_outfile,data))
		print data,'....finished!!!'
###########################################################################################################################3
if __name__=="__main__":
	bedfile_transform()
	run_VEP_1000genome(path_soft_VEP,path_VEP)	
	from_VEP_to_FATHMM_chunk()
	fathmm_process_VEP(path_soft_fathmm)
