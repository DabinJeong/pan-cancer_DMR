#!/usr/local/bin/python3

#import pandas as pd
#import numpy as np
#
##Parse ppi_nwk by genes that are in mart_genes
#ppi_nwk = pd.read_csv("9606.protein.links.full.v10.5.txt", sep=' ',header=0,index_col=None)
#mart = pd.read_csv("mart_export.txt", sep='\t', header=0, index_col=None)
#
#stripped_ppi_column1 = [string.lstrip("9606.") for string in ppi_nwk['protein1']]
#stripped_ppi_column2 = [string.lstrip("9606.") for string in ppi_nwk['protein2']]
#ppi_nwk[['protein1']] = stripped_ppi_column1
#ppi_nwk[['protein2']] = stripped_ppi_column2
#
#index = np.logical_and(ppi_nwk.iloc[:,0].isin(mart.iloc[:,3]), ppi_nwk.iloc[:,1].isin(mart.iloc[:,3]))
#
#filtered_nwk = ppi_nwk[index]
#
#print( ppi_nwk.shape )
#print(filtered_nwk.shape)

import pandas as pd
import numpy as np

'''Parse ppi_nwk by genes that are in mart_genes'''

ppi_nwk = pd.read_csv("9606.protein.links.full.v10.5.txt", sep=' ',header=0,index_col=None)
mart = pd.read_csv("mart_export.txt", sep='\t', header=0, index_col=None)

#Cleaning data
stripped_ppi_column1 = [string.lstrip("9606.") for string in ppi_nwk['protein1']]
stripped_ppi_column2 = [string.lstrip("9606.") for string in ppi_nwk['protein2']]
ppi_nwk[['protein1']] = stripped_ppi_column1
ppi_nwk[['protein2']] = stripped_ppi_column2

#Parse ppi_nwk by mart genes
index = np.logical_and(ppi_nwk.iloc[:,0].isin(mart.iloc[:,3]), ppi_nwk.iloc[:,1].isin(mart.iloc[:,3]))
filtered_nwk = ppi_nwk[index]

#Substitute proteinIDs with gene name
# building gene name distionary (keys: protien ID, values: gene symbol)
keys = [str(i) for i in mart['Protein stable ID']]
values = [str(f) for f in mart['Gene name']]
dictionary = dict(zip(keys, values))
dictionary.pop('nan',None)
# apply the conversion to the filtered network
filtered_nwk.iloc[:,[0]] = [dictionary[i] for i in filtered_nwk['protein1']]
filtered_nwk.iloc[:,[1]] = [dictionary[i] for i in filtered_nwk['protein2']]

#Write the network to a file
filtered_nwk.to_csv('filtered_ppi.nwk', sep='\t', index=None, header=True)
