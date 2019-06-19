#!/usr/bin/env python
import sys
import argparse
import scipy.stats as stats
import pandas as pd

parser=argparse.ArgumentParser(
        usage='''\
%(prog)s [options] gene2subtype -trait2genes trait2genes 

example: %(prog)s gene2subtype.txt -sym2entrez sym2entrez.txt -trait2genes trait2genes.txt -backgroundGenes backgroundGemes.txt -pcut 0.05 -topK None -o outfile.txt
''')

parser.add_argument('gene2subtype', metavar='str', help='gene2subtype file')
parser.add_argument('-sym2entrez',required=False, help='symbol to entrez ID mapped file')
parser.add_argument('-trait2genes', required=False, metavar='str', default='/data1/project/hongryul/GOanalysis/GOBPname2gene.arabidopsis.BP.20180907.txt', help='trait2genes file')
parser.add_argument('-backgroundGenes', required=False, metavar='str', default='None', help='allgenes in first column file')
parser.add_argument('-method', required=False, metavar='[fisher|binomial]', default='fisher', help='method for statistical test')
parser.add_argument('-s', required=False, type=int, metavar='N', default=0, help='skipped lines for allgenes file')
parser.add_argument('-pcut', required=False, type=float, metavar='N', default=1.0, help='pvalue cutoff')
parser.add_argument('-topK', required=False, type=int, metavar='N', default=None, help='show top K result')
parser.add_argument('-o', dest='outfile', required=False, metavar='str', default='stdout', help='outfile')
args=parser.parse_args()

if args.outfile == 'stdout':
	OF=sys.stdout
else:
	OF=open(args.outfile,'w')

if args.sym2entrez != 'None':
	sym2entrez=pd.read_csv(args.sym2entrez,sep='\t',header=None,index_col=0).to_dict()[1]
	
if args.backgroundGenes != 'None':
	set_gene=set()
	IF=open(args.backgroundGenes,'r')
	for line in IF:
		gene=line.rstrip('\n').split('\t',1)[0]
		set_gene.add(gene)

dic_gene2trait={}
dic_trait2gene={}
IF=open(args.trait2genes,'r')
for line in IF:
	trait,genes=line.rstrip('\n').split('\t')
	lst_gene=genes.split(',')
	for gene in lst_gene:
		if args.backgroundGenes != 'None' and gene not in set_gene:
			continue
		if trait not in dic_trait2gene:
			dic_trait2gene[trait]=set()
		if gene not in dic_gene2trait:
			dic_gene2trait[gene]=set()
		dic_gene2trait[gene].add(trait)
		dic_trait2gene[trait].add(gene)

dic_trait2ratio={}
for trait in dic_trait2gene.keys():
	dic_trait2ratio[trait]=float(len(dic_trait2gene[trait]))/len(dic_gene2trait)

	
def GO_enrichment(lst_gene):
	dic_trait2count={}
	set_tested_gene=set()
	for gid in lst_gene:
		if gid not in dic_gene2trait:
			continue
		for trait in dic_gene2trait[gid]:
			if not trait in dic_trait2count:
				dic_trait2count[trait]=0
			dic_trait2count[trait]+=1
		set_tested_gene.add(gid)
	lst_out=[]
	for trait, count in dic_trait2count.items():
		occured_in_tested, total_tested, occured_in_background, total_background = count, len(set_tested_gene), len(dic_trait2gene[trait]), len(dic_gene2trait)
		if occured_in_tested == 0:
			pval=1.0
		else:
			if args.method == 'binomial':
				pval=1.0-stats.binom.cdf(occured_in_tested-1, total_tested, dic_trait2ratio[trait])	# ovccured_in_test-1 means p(X>=n) i.e. contain
			elif args.method == 'fisher':
				oddratio,pval=stats.fisher_exact([[occured_in_tested, total_tested-occured_in_tested], [occured_in_background-occured_in_tested, total_background-total_tested-occured_in_background+occured_in_tested]], alternative='greater')	# 2X2 fisher's exact test
		lst_out.append([trait, pval, occured_in_tested, total_tested, occured_in_background, total_background])
	return sorted(lst_out, key=lambda x: x[1])

IF=open(args.gene2subtype,'r')
for i in range(args.s):
	IF.readline()
OF.write('\t'.join(['subtype','setid','pval','#occured_in_tested','#total_tested','#occured_in_background','#total_background'])+'\n')
dic_subtype2gene={}
for line in IF:
	s=line.rstrip().split('\t')
	if args.sym2entrez == 'None':
		if len(s) == 1:
			gene,subtype=s[0],'None'
		else:
			gene,subtype=s[0:2]
		if subtype not in dic_subtype2gene:
			dic_subtype2gene[subtype]=set()
		dic_subtype2gene[subtype].add(gene)
	else:
		if len(s) == 1:
			gene,subtype=str(sym2entrez[s[0]]),'None'
		else:
			gene,subtype=str(sym2entrez[s[0]]), s[1]
		if subtype not in dic_subtype2gene:
			dic_subtype2gene[subtype]=set()
		dic_subtype2gene[subtype].add(gene)
for subtype, lst_gene in sorted(dic_subtype2gene.items(),key=lambda x:float(x[0]) if x[0].isdigit() else x[0]):
	lst_out=GO_enrichment(lst_gene)
	if len(lst_out) == 0:
		OF.write('\t'.join(map(str, [subtype]+['nan']*6))+'\n')
	if args.topK == None:
		for out in lst_out:
			if float(out[1]) <= args.pcut:
				OF.write('\t'.join(map(str, [subtype, out[0]]+out[1:]))+'\n')
	else:
		for out in lst_out[0:args.topK]:
			if float(out[1]) <= args.pcut:
				OF.write('\t'.join(map(str, [subtype, out[0]]+out[1:]))+'\n')
