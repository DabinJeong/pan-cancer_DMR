import pandas as pd
import argparse
import scipy.stats as stats
from math import log10,log10
import pandas as pd
import sys

parser=argparse.ArgumentParser(usage='')
parser.add_argument('-methyl',required=True, help='')
parser.add_argument('-community',required=True, help='')
#parser.add_argument('-total_node',required=True, help='')
parser.add_argument('-o',help='')
args=parser.parse_args()

df_methyl_0 = pd.read_csv(args.methyl,sep='\t')
df_comm = pd.read_csv(args.community,sep='\t')
total_node = 18208 

df_methyl = df_comm.merge(df_methyl_0,on=['gene'],how="inner").rename(columns={'value':'methylation_value'}).groupby(['gene','community']).mean().reset_index()

df_methyl['methyl_flag'] = df_methyl['methylation_value']>0
df_methyl_cnt = df_methyl.groupby(['community','methyl_flag']).size().unstack(fill_value=0)
for flag in [True,False]:
    if flag not in df_methyl_cnt.columns.tolist():
        df_methyl_cnt[flag]=0

df_methyl_cnt.rename(columns={False:'methyl-',True:'methyl+'},inplace=True)
L=df_methyl.shape[0]
#L_plus, L_minus = df_methyl_cnt['methyl+'].sum(), df_methyl_cnt['methyl-'].sum()
df_methyl_cnt['#DMR']=df_methyl_cnt['methyl+']+df_methyl_cnt['methyl-']
df = df_methyl_cnt.join(df_comm.groupby(['community']).count().rename(columns={'gene':'gene_cnt'}))

#a_plus=df[True], l=df['gene_cnt']
#df['methyl+pval']=df.apply(lambda x: stats.fisher_exact([[x['methyl+'],L_plus-x['methyl+']],[x.gene_cnt-x['methyl+'],total_node-(L_plus+x.gene_cnt-x['methyl+'])]])[1],axis=1)
#df['methyl-pval']=df.apply(lambda x: stats.fisher_exact([[x['methyl-'],L_plus-x['methyl-']],[x.gene_cnt-x['methyl-'],total_node-(L_plus+x.gene_cnt-x['methyl-'])]])[1],axis=1)
df['DMRpval']=df.apply(lambda x: stats.fisher_exact([[x['#DMR'],L-x['#DMR']],[x.gene_cnt-x['#DMR'],total_node-(L+x.gene_cnt-x['#DMR'])]])[1],axis=1)
df['-log(DMRpval)']=df['DMRpval'].apply(lambda x:-log10(x) if x>0 else 1)
#df['-log(methyl+/methyl-)']=df['-log(methyl+/methyl-)'].apply(lambda x:-log10(x))
#df.reset_index()[['community','methyl+pval','methyl-pval']].to_csv(args.o, sep='\t',index=False)
df[df['DMRpval']<0.05].reset_index()[['community','methyl+','methyl-','-log(DMRpval)','DMRpval']].to_csv(args.o, sep='\t',index=False)
