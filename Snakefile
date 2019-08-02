cancerLi = "BRCA LAML".split()

nwkCorr_corrCut="0.5"
nwkCorr_FCcut="0.15"
GO_pCut="0.05"
ppi="String_db/filtered_ppi.nwk"

rule all:
    input:
        "jaccard_result"

rule nwk_corr:
    input:  
        ppi={ppi},
        exp="data_processed/{cancer}.mRNA.TPM.selected.onlyCancer.txt",
        log2FC="data_processed/log2FC_files/{cancer}_log2FC.txt",
        mutSample="data_processed/{cancer}.mutated_sampleID",
    output:
        "{cancer}.filtered.noExpCut.sorted.nwk"        
    log:
        "log.{cancer}.nwk_corr"
    shell:
        """
        python assign_score.py {input.ppi} {input.exp} {input.log2FC} -mutSamples {input.mutSample} -cancerType {wildcards.cancer} -corrCut {nwkCorr_corrCut} -FCcut {nwkCorr_FCcut} -o {output}
        """

rule nwk_module:
    input:
        nwk="{cancer}.filtered.noExpCut.sorted.nwk"
    output:
        comm="{cancer}.community",
        edge="{cancer}.community.filtered_edgelist"
    log:
        "log.{cancer}.nwk_module"
    shell:
        """
        Rscript community_detection_mlc_v3.R {input.nwk} {output.comm} {output.edge}
        """

rule GO:
    input:
        comm="{cancer}.community"
    output:
       GO="{cancer}.GO" 
    log:
       "log.{cancer}.GO"
    shell:
        """
        python ./myGOenrichment_DB.py {input.comm} -sym2entrez GO_score/sym2entrez/{wildcards.cancer}.* -trait2genes GO_score/GOBPname2genes.human.BP.txt -pcut {GO_pCut} -s 1 -topK 1 -o {output.GO} 
        """

rule ttest:
    input:
        comm="{cancer}.community",
        log2FC="data_processed/log2FC_files/{cancer}_log2FC.txt",
        edge_filtered="{cancer}.community.filtered_edgelist",
        GO="{cancer}.GO"
    output:
        comm2="{cancer}.ttest_go_pval_cut_community",
        edge2="{cancer}.ttest_go_pval_cut_edgelist"
    log:
        "log.{cancer}.ttest"
    shell:
        """
        python 1sample_ttest.py -c {input.comm} -l {input.log2FC} -nwk {input.edge_filtered} -go {input.GO} -cancer {wildcards.cancer} -o1 {output.comm2} -o2 {output.edge2} >& {log}
        """

rule methyl_promoter:
    input:
        methyl="methylation/{cancer}_DMR.txt"
    output:
        mehtyl_P="methylation/{cancer}_promoterDMR.txt"     
    run:
        import pandas as pd 
        df=pd.read_csv(input.methyl,sep='\t') 
        df[df['region']=='promoter'][['name','value']].rename(columns={'name':'gene'}).to_csv(output.mehtyl_P,index=False,sep='\t')

rule fisher:
    input:
        methyl="methylation/{cancer}_promoterDMR.txt",
        comm="{cancer}.ttest_go_pval_cut_community"
    output:
        "{cancer}.fisher.out"
    shell:
        """
        python fisher.py -methyl {input.methyl} -community {input.comm} -o {output}
        """

rule cnt_gene:
    input:
        cancer_comm = "{cancer}.ttest_go_pval_cut_community"
    params:
        cancer = "{cancer}"
    output:
        cancer_cnt = temp("{cancer}.gene_cnt_per_cluster")
    run:
        import pandas as pd 
        d={}
        df=pd.read_csv(input.cancer_comm,sep='\t')
        df2=pd.DataFrame(df.groupby('community')['gene'].apply(list))
        for comID in df2.index:
            d['{}_{}'.format(params.cancer,comID)] = [df2.loc[comID,'gene']]
        pd.DataFrame.from_dict(d, orient='index',columns=['geneList']).to_csv(output.cancer_cnt,header=False,sep='\t')    

rule concat_cnt:
    input:
        expand("{cancer}.gene_cnt_per_cluster", cancer=cancerLi)
    output:
        temp("gene_cnt_per_cluster")
    shell:
        """
        cat {input} >> {output}
        """

rule jaccard:
    input:
        pan_cnt = "gene_cnt_per_cluster"
    output:
        jac = "jaccard_result"
    run:
        import pandas as pd
        import numpy as np
        from itertools import combinations
        d1=pd.read_csv(input.pan_cnt,sep='\t',header=None).set_index(0).T.to_dict('list')
        d = {k:v[0].strip("[]").split(',') for k,v in d1.items()}
        print(d)
        nC2=combinations(d.keys(),2)
        with open(output.jac,'w') as f:
            for pair in nC2:
                set1,set2=d[pair[0]],d[pair[1]]
                j_idx=float(np.intersect1d(set1,set2).size)/float(np.union1d(set1,set2).size)
                print ("{}   {}  {}".format(pair[0],pair[1],j_idx), file=f)

rule clear:
    input:
        expand("{cancer}.GO",cancer=cancerLi)
    shell:
        """rm -f *.community* |rm *.GO|rm *.filtered.noExpCut.*| rm -f log.*|rm -f *.ttest_go_pval_cut_community|rm -f *.ttest_go_pval_cut_edgelist"""
