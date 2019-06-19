import pandas as pd
import argparse

parser=argparse.ArgumentParser(usage='%(prog)s sym2entrez_dict community_profile out')
parser.add_argument('sym2entrez_dict')
parser.add_argument('gene_list')
parser.add_argument('out')
args=parser.parse_args()
	
sym2entrez=pd.read_csv(args.sym2entrez_dict,sep='\t',header=None,index_col=0).to_dict()[1]
df=pd.read_csv(args.gene_list, sep='\t',names=['gene'])
df['gene']=df['gene'].apply(lambda x:sym2entrez[x])
df.to_csv(args.out,sep='\t',index=False,header=False)
