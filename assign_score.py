import pandas as pd
import argparse
import numpy as np
import time
import sys
import multiprocessing
from scipy.stats.stats import pearsonr

def isNum(x):
    try:
        float(x)
        return True
    except:
        return False

def corrCut(nwk,cutoff=None):
    '''correlation cutoff, positive sorting'''
    nwk.dropna(subset=['string*corr'],inplace=True)
    nwk.sort_values(by=['string*corr'],inplace=True,ascending=False)
    if cutoff!=None: 
        return nwk[nwk['correlation']>=cutoff][['protein1','protein2','string*corr']]
    else:
        return nwk[['protein1','protein2','string*corr']]

def setMinExp(nwk,exp,expCut):
    '''remove gene whose expression is lower than expCut'''
    ###nwk = pd.read_csv(tableFile, sep='\t',header=0,index_col=None)
    filtered_gene = exp[exp.max(axis=1)>1]['Hybridization REF']
    boolean_mask = np.logical_and(nwk['protein1'].isin(filtered_gene),nwk['protein2'].isin(filtered_gene))
    return nwk[boolean_mask]

def expCut(nwk,exp,sample_list,expCut):
    '''remove gene whose mean of group(mutated/not-mutated) expression is lower than expCut'''
    with open(sample_list) as f:
        mutSamples=f.read().strip().split()
    exp['no_mut']=exp.loc[:,~exp.columns.isin(mutSamples)].mean(axis=1)
    exp['mut']=exp.loc[:,exp.columns.isin(mutSamples)].mean(axis=1)
    boolean_mask = np.logical_or(exp['no_mut']>=1,exp['mut']>=1)
    gene_selected = exp[boolean_mask]['Hybridization REF']
    boolean_mask2 = np.logical_and(nwk['protein1'].isin(gene_selected),nwk['protein2'].isin(gene_selected))
    return nwk[boolean_mask2]

def FCcut(nwk,FC_df,FCcut):
    ###nwk = pd.read_csv(tableFile, sep='\t',header=0,index_col=None)
    #building dictionary
    keys = FC_df.iloc[:,0]
    values = FC_df.iloc[:,1]
    dictionary = dict(zip(keys, values))
    dictionary.pop('?','Not Found')
    first_col = np.array([dictionary[i] for i in nwk.iloc[:,0]])
    second_col = np.array([dictionary[i] for i in nwk.iloc[:,1]])
    boolean_mask = np.logical_and(abs(first_col) >= FCcut, abs(second_col) >= FCcut)
    boolean_mask2 = nwk['protein1'].apply(lambda x: dictionary[x])*nwk['protein2'].apply(lambda x: dictionary[x])>0
    bigFC_nwk = nwk[boolean_mask & boolean_mask2]
    #bigFC_nwk.to_csv(outFile, sep='\t', index=None, header=True)
    return bigFC_nwk

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='python %(prog)s nwkFile seedFile -cancerType [BRCA] -FC [log2FC.txt]')
    parser.add_argument('STRING',help='STRING')
    parser.add_argument('exp',help='exp File')
    parser.add_argument('log2FC',help='log2FC file $cancerType_log2FC.txt')
    parser.add_argument('-mutSamples',required=False,help='DNA methyl transferse mutated sample list (TCGA ID)')
    parser.add_argument('-cancerType',required=True,help='TCGA cancer ID')
    parser.add_argument('-corrCut',type=float, required=False,help='corealtion cutoff')
    parser.add_argument('-FCcut',type=float, required='log2FC' in sys.argv, help='log2FC cutoff')
    #parser.add_argument('-expCut',type=float, required='exp' in sys.argv, help='exp mean cutoff')
    parser.add_argument('-o',required=True,help='output')
    args=parser.parse_args()
        
    ####correaltion score combined with string score
    start=time.time()

    exp=pd.read_csv(args.exp,sep='\t',header=0)
    
    log2FC=pd.read_csv(args.log2FC, sep='\t', header=0, index_col=None)
    #remove duplicates
    data=[]
    with open(args.STRING) as f:
        for line in f.readlines():
            tmp=line.strip().split()
            data.append(sorted(tmp[:2])+tmp[2:])
    STRING=pd.DataFrame(data[1:],columns=data[0],dtype=np.float)
    STRING.drop_duplicates(subset=['protein1','protein2'],inplace=True)
    mask=np.logical_and(STRING['protein1'].isin(exp['Hybridization REF']),STRING['protein2'].isin(exp['Hybridization REF']))
    STRING=STRING[mask]
    #filter network with expCut, FCcut
    #STRING_expCut=expCut(STRING,exp,args.mutSamples,args.expCut)
    STRING_FCcut=FCcut(STRING,log2FC,args.FCcut)
    #STRING_expCut_FCcut=FCcut(STRING_expCut,log2FC,args.FCcut)
    
    #make exp dictionary to calculate correlation
    lst_exps=dict() 
    with open(args.exp) as f:
        lines=f.readlines()
    for line in lines:
        s=line.strip().split('\t')
        if not isNum(s[1]):
            continue
        else:
            gene, exps = s[0], list(map(float,s[1:]))
            lst_exps[gene]=exps
    #lst_pairs=zip(STRING_expCut_FCcut['protein1'],STRING_expCut_FCcut['protein2'])
    lst_pairs2=zip(STRING_FCcut['protein1'],STRING_FCcut['protein2'])
    def myCorr(x):
        g1,g2=sorted(x)
        if g1 == g2:
            val = 0.0
        else:
            val, pval = pearsonr(lst_exps[g1],lst_exps[g2])
        return (g1,g2,val)

    pool = multiprocessing.Pool(15)    
    #res=pool.imap_unordered(myCorr, lst_pairs)
    res2=pool.imap_unordered(myCorr, lst_pairs2)
    #corr_res=[]
    corr_res2=[]
    #for g1,g2,val in res:
    #    if g1==g2:
    #        continue
    #    corr_res.append([g1,g2,val])
    for g1,g2,val in res2:
        if g1==g2:
            continue
        corr_res2.append([g1,g2,val])
    #corr_dict=pd.DataFrame(corr_res,columns=['protein1','protein2','correlation'])
    corr_dict2=pd.DataFrame(corr_res2,columns=['protein1','protein2','correlation'])
    
    #STRING_expCut_FCcut=pd.merge(corr_dict,STRING_expCut_FCcut,on=['protein1','protein2'])
    #STRING_expCut_FCcut['combined_score_norm']=STRING_expCut_FCcut['combined_score'].apply(lambda x:float(x)/1000)
    #STRING_expCut_FCcut['string*corr']=STRING_expCut_FCcut['correlation']*STRING_expCut_FCcut['combined_score_norm']
    #STRING_expCut_FCcut_corrCut=corrCut(STRING_expCut_FCcut,args.corrCut)
    #STRING_expCut_FCcut_corrCut[['protein1','protein2','string*corr']].to_csv(args.cancerType+'.filtered.sorted',index=False,sep='\t') 
    
    STRING_FCcut=pd.merge(corr_dict2,STRING_FCcut,on=['protein1','protein2'])
    STRING_FCcut['combined_score_norm']=STRING_FCcut['combined_score'].apply(lambda x:float(x)/1000)
    STRING_FCcut['string*corr']=STRING_FCcut['correlation']*STRING_FCcut['combined_score_norm']
    STRING_FCcut_corrCut=corrCut(STRING_FCcut,args.corrCut)
    #STRING_FCcut_corrCut[['protein1','protein2','string*corr']].to_csv(args.cancerType+'.filtered.noExpCut.sorted.DNMT3A_only',index=False,sep='\t')
    STRING_FCcut_corrCut[['protein1','protein2','string*corr']].to_csv(args.o, index=False,sep='\t')
    end=time.time()
    time_elapsed=end-start
    print ('network_correaltion: '+str(time.strftime("%H:%M:%S",time.gmtime(time_elapsed))))
   

    ####filter data
    #df=corrCut(args.cancerType+'.table',args.corrCut)
    #df2=FCcut(df,args.log2FC,args.FCcut)
    #df3=setMinExp(df2,args.exp,1)
    #df3=expCut(df,args.exp,args.mutSamples,1)
    #df3.to_csv(args.cancerType+'.corr.'+ str(args.corrCut) +'.FC.'+str(args.FCcut)+'.expCut',sep='\t',index=False)
