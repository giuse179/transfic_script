#!/bin/python
#12/04/2013

import sys, csv, re
import pandas as pd
import numpy as np
#############################################################################################
path_fathmm=all_fathhmm_sorted
path_go=prova_50000'		
#########################################################################################################################
fileout1=open('baseline_tollerance','w')
##########################################################################
def parse_fathmm(path):					
	dataframes=pd.read_table(path,index_col=2,na_values=['No Prediction Available','No Sequence Record Found','Score'])	
	cols=['Score','Substitution']		
	dataframes=dataframes.ix[:,cols].copy()									
	return dataframes

def parse_go(path):
	cols=['Gene_ID','Protein_ID','GO_Accession']	
	all_data=pd.read_csv(path,header=None,index_col=1,skiprows=1,names=cols,na_values='').drop_duplicates()
	grouped=all_data.groupby(['GO_Accession'])
	return grouped,all_data

def global_mean(fathmm,all_go):	
	data=fathmm.join(all_go,how='inner').sort_index(by=['GO_Accession']).drop_duplicates()
	all_mean=data['Score'].mean()
	all_std=data['Score'].std()
	data=pd.DataFrame()
	return all_mean,all_std 
		
def merge(fathmm,go):
	b=fathmm
	a=go	
	c=b.join(a,how='inner').sort_index(by=['GO_Accession'])	
	return c
		
def statistic_goGroup(data):
	go_group=data.groupby(['GO_Accession'])
	mean=data.join(go_group['Score'].mean(),on='GO_Accession',rsuffix='_mean',how='inner')
	all_=mean.join(go_group['Score'].std(),on='GO_Accession',rsuffix='_std',how='inner')	
	return all_

def count_(data,all_mean,all_std):
	data['mean_TOT']=all_mean
	data['std_TOT']=all_std
	data=data.reset_index()
	data['Protein_ID']=data['index']
	col=['GO_Accession','Gene_ID','Substitution','Protein_ID','Score_mean','Score_std','mean_TOT','std_TOT']	
	data=data.ix[:,col].copy()
	return data.drop_duplicates()

def groped(data):
	Gene_go_group=data.groupby(['GO_Accession','Gene_ID','Substitution'])
	agg_count=pd.DataFrame(Gene_go_group.size().fillna(np.nan))
	agg_count=agg_count.reset_index()
	if not agg_count.empty:
		matrix_info= pd.DataFrame(agg_count['Gene_ID'].value_counts(),columns=['variants_for_Gene_count'])
		matrix_info['variants_count']=agg_count['GO_Accession'].value_counts()[0]
		matrix_info['gene_count']=matrix_info.count()[0]
		matrix_info['info']= np.where(matrix_info['variants_count'] < 20,'ALL','POOLED')
		matrix_info['Gene_ID']=matrix_info.index
		return matrix_info
	
def table(data,info):
	if info is not None:	
		newdata=pd.merge(data,info)
		return newdata

def final_mean(data):
	if data is not None:
		data['final_mean']=np.where(data['info']=='ALL',data['mean_TOT'],data['Score_mean'])
		data['final_std']= np.where(data['info']=='ALL',data['std_TOT'],data['Score_std'])
		cols=['GO_Accession','Gene_ID','Protein_ID','final_mean','final_std','variants_count','info']
		z=data.ix[:,cols].sort_index(by=['Gene_ID'])		
		z=z.drop_duplicates().copy()
		return z

def concatenated(data,all_data):		
	if data is not None:
		all_data.append(data[1:])
	df=pd.concat(all_data)
	return df

def matrix(data):
	data=data.sort_index(by=['Gene_ID','variants_count'])
	data=data.drop(['GO_Accession'],axis=1)
	group=pd.DataFrame(data.groupby(['Gene_ID']).sort_index(by=['Gene_ID','variants_count']))
	grouped=pd.DataFrame(data.groupby(['Gene_ID']).size(),columns=['Gene_multiple'])
	grouped2=pd.DataFrame(data.groupby(['Gene_ID','info']).size(),columns=['Gene_Info_multiple'])
	grouped3=pd.DataFrame(data.groupby(['Gene_ID','info'])['variants_count'].min(),columns=['min_value'],dtype='int64')
	data1=group.reset_index()
	data2=grouped.reset_index()
	data3=grouped2.reset_index()
	data3_b=grouped3.reset_index()
	c=pd.merge(data1,data2,on='Gene_ID')
	data4=pd.merge(c,data3,on=['Gene_ID','info'])
	data5=pd.merge(data4,data3_b,on=['Gene_ID','info'])
	data5=data5.drop(['index'],axis=1)
	data5=data5.set_index('Gene_ID').drop_duplicates()
	return data5

def unique(data):
	data1=data[(data['variants_count']==data['min_value'])]   # & (data['info']=='POOLED')]
	data1['unique']=np.where((data1['Gene_multiple']>1)&(data1['info']=='POOLED'),1,np.where(data1['Gene_multiple']==1,1,np.where(data1['Gene_multiple']==data1['Gene_Info_multiple'],1,0)))
	data2=data1[data1['unique']==1]
	data2=data2.drop(['variants_count','Gene_multiple','Gene_Info_multiple','min_value','unique'],axis=1)
	return data2


############### RUN ###########################
if __name__=="__main__":					
	all_=[]	
	i=0
	FATH_HMM= parse_fathmm(path_fathmm)
	grouped_GO,all_data= parse_go(path_go)
	mean, std=global_mean(FATH_HMM,all_data)
	print 'start'
	for name,group in grouped_GO:		
		i+=1	
		data=count_(statistic_goGroup(merge(FATH_HMM,group)),mean,std)
		info=groped(data)
		r=table(data,info)
		final=final_mean(r)		
		alls=concatenated(final,all_)
		print i
	p=matrix(alls)
	print 'matrix'
	final_all=unique(p)
	final_all.to_csv(fileout1,sep='\t',index=True,header=True) #,cols=['GO_Accession','Gene_ID_x','Protein_ID','final_mean','final_std','info']	






















