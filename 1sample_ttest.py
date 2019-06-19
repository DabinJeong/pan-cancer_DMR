#!/usr/local/bin/python3

'''performing 1sample ttest for each community of a cancer network'''

##Imports
import pandas as pd
import numpy as np
from scipy import stats
import argparse

# Argument parsing
parser = argparse.ArgumentParser(description='')
parser.add_argument('-c','--community', type=str, metavar='', required=True, help='community_info ${cluster}/$[cancer].corr.0.5_FC.0.25_community')
parser.add_argument('-l', '--logFC', type=str, metavar='', required=True, help='log2FC file $cancerType_log2FC.txt')
parser.add_argument('-nwk', type=str, metavar='', required=True, help='edgelist for parsing')
parser.add_argument('-go',type=str,required=True)
parser.add_argument('-cancer', type=str, metavar='', required=True, help='cancer type')
parser.add_argument('-o1',type=str,required=True,help="community output file name")
parser.add_argument('-o2',type=str,required=True,help="edgelist output file name")
args = parser.parse_args()

# Import Data
clusters = pd.read_csv(args.community, sep='\t',header=0,index_col=None)
FC_df = pd.read_csv(args.logFC, sep='\t', header=0, index_col=None)
GO_df = pd.read_csv(args.go,sep='\t')

keys = FC_df.iloc[:,0]
values = FC_df.iloc[:,1]
dictionary = dict(zip(keys, values))
dictionary.pop('?','Not Found')

FC_col = [ dictionary[i] for i in clusters.iloc[:,0] ]

new_data=pd.DataFrame(np.column_stack([clusters,FC_col]))

new_data.columns = ['gene', 'community', 'log2FC']

# finding commmunity with only 1 member
# exclude community with only 1 member from the ttest input arrays
community_numbers = new_data.community.unique()[new_data.groupby('community').count().gene>=2]

## ttest_results = [['community : {}'.format(i), stats.ttest_1samp(new_data[new_data.community==i].log2FC , 0, 0)[0],
##     stats.ttest_1samp(new_data[new_data.community==i].log2FC , 0, 0)[1]]
##     for i in community_numbers]

ttest_results = [[i, stats.ttest_1samp(new_data[new_data.community==i].log2FC , 0, 0)[0],
    stats.ttest_1samp(new_data[new_data.community==i].log2FC , 0, 0)[1]]
    for i in community_numbers]

concat_result = pd.DataFrame(ttest_results)
concat_result.columns = ['community_num', 'test_result', 'p-value']
concat_result = concat_result.sort_values(by=['p-value'])

# Fetch only the p-values smaller than 10^(-8) (+GO p-value)
concat_result_sorted = concat_result[concat_result['p-value']<=10**(-9)]
community_ttest = set(concat_result[concat_result['p-value']<=10**(-9)]['community_num'].tolist())
community_go = set(GO_df[GO_df['pval']<10**(-9)]['subtype'].tolist())
community_inter = community_ttest.intersection(community_go)
new_data[new_data['community'].isin(community_inter)][['gene','community']].to_csv(args.o1,sep='\t',index=None, header=True)
#new_data[new_data['community'].isin(community_inter)][['gene','community']].to_csv(args.cancer+'.ttest_go_pval_cut_community.DNMT3A_only',sep='\t',index=None, header=True)

# Prepare edgelist for parsing (by small p-valued communities)
edge_list = pd.read_csv(args.nwk, sep='\t',header=0,index_col=None)

# Dictionary : (key : gene name. value : community_num)
cluster_keys = clusters.iloc[:,0]
cluster_values = clusters.iloc[:,1]
cluster_dict = dict(zip(cluster_keys, cluster_values))

clusters_source = pd.Series([cluster_dict[i] for i in edge_list.iloc[:,0]])
clusters_target = pd.Series([cluster_dict[i] for i in edge_list.iloc[:,1]])

# Check if both source and the target node is in the selected communities
# Parse Edge_list by the selected communities
#bool_mask = np.logical_and(clusters_source.isin(concat_result_sorted['community_num']), clusters_target.isin(concat_result_sorted['community_num']))
bool_mask = np.logical_and(clusters_source.isin(community_inter),clusters_target.isin(community_inter))
edge_list_parsed = edge_list[bool_mask]

#edge_list_parsed_filename = ''.join([args.cancer,'.ttest_go_pval_cut_edgelist.DNMT3A_only'])
edge_list_parsed.to_csv(args.o2, sep='\t', index=None, header=True)
