import pandas as pd
import argparse
from math import log10
import numpy as np
import sys
import time

def funcScore(df_Comm,df_GO):
	'''
	Calculate functionality score for each community	
	================================================
	Input
		df_Comm (pd.DataFrame)
		df_GO(pd.DataFrame)
	Output
		funcScore (flaot)
		'''
	sizeComm=df_Comm.groupby(['community']).size().reset_index(name='size')
	sizeTotal=df_Comm.shape[0]
	sizeComm['weight']=sizeComm['size']/sizeTotal
	cntComm=sizeComm.shape[0]
	df_GO['pval'].fillna(1,inplace=True)
	df_GO['pval'].replace(0,10**(-300),inplace=True)
	sizeComm['funcScore']=-df_GO['pval'].apply(lambda x:log10(x))*sizeComm['weight']
	return sizeComm['funcScore'].sum()


if __name__=="__main__":
	parser=argparse.ArgumentParser(usage='%(prog)s communityFile GOfile -topN 1000')
	parser.add_argument('communityFile')
	parser.add_argument('GOfile')
	parser.add_argument('-topN',required=True)
#	parser.add_argument('-o')
	args=parser.parse_args()
	df1=pd.read_csv(args.communityFile,sep='\t')
	df2=pd.read_csv(args.GOfile,sep='\t')
	print '{}\t{}'.format(args.topN, funcScore(df1,df2))
