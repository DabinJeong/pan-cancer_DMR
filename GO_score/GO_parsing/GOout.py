#!/usr/bin/env python
import sys
import argparse

parser=argparse.ArgumentParser(
	usage='''\
%(prog)s [options] gene2GO go.obo

example: %(prog)s gene2GO go.obo
''')

parser.add_argument('gene2GO', help='gene2GO file (1st column:gene, 2nd column:GO')
parser.add_argument('obo', help='go.obo file')
parser.add_argument('-namespace', required=False, default='BP', help='BP: biological process, CC: cellular component, MF: molecular function, ALL: all')
parser.add_argument('-outfmt', required=False, default='GOname2genes', help='gene2GO, GOname2genes, etc.')
parser.add_argument('-o', dest='outfile', required=False, metavar='str', default='stdout', help='outfile')
args=parser.parse_args()

if args.outfile == 'stdout':
	OF=sys.stdout
else:
	OF=open(args.outfile,'w')

dic_GO2name={}
dic_GO2namespace={}
dic_alt_id2GO={}
dic_term2content={}
IF=open(args.obo,'r')
for line in IF:
	if line.startswith('[Term]') or line.startswith('[Typedef]'):
		if line.startswith('[Term]') and 'id' in dic_term2content:
			if len(dic_term2content['id']) != 1: print 'obo parse error: # id is over 1'; sys.exit()
			id=dic_term2content['id'][0]

			if len(dic_term2content['name']) != 1: print 'obo parse error: # name is over 1'; sys.exit()
			name=dic_term2content['name'][0]

			if len(dic_term2content['namespace']) != 1: print 'obo parse error: # namespace is over 1'; sys.exit()
			namespace=dic_term2content['namespace'][0]
			if namespace == 'biological_process': namespace = 'BP'
			elif namespace == 'cellular_component': namespace = 'CC'
			elif namespace == 'molecular_function': namespace = 'MF'
			else:  print 'obo parse error: wrong namespace' + namespace; sys.exit()

			if 'alt_id' in dic_term2content:
				for alt_id in dic_term2content['alt_id']:
					dic_alt_id2GO[alt_id] = id
			dic_GO2name[id]=name
			dic_GO2namespace[id]=namespace
		dic_term2content={}
		continue
	s=line.rstrip().split(':',1)
	if len(s) < 2:
		continue
	term, content= s[0].strip(), s[1].strip()
	if term not in dic_term2content:
		dic_term2content[term]=[]
	dic_term2content[term].append(content)

dic_gene2GOs={}
dic_GO2genes={}
IF=open(args.gene2GO,'r')
for line in IF:
	if line.startswith('#'):
		continue
	s=line.rstrip('\n').split('\t')
	gene, GO = s[0], s[1]
	if GO in dic_alt_id2GO:
		GO = dic_alt_id2GO[GO]
	if GO not in dic_GO2namespace:
		continue
	namespace = dic_GO2namespace[GO]
	if args.namespace != 'ALL' and namespace != args.namespace:
		continue

	if GO not in dic_GO2genes:
		dic_GO2genes[GO] = set()
	dic_GO2genes[GO].add(gene)

	if gene not in dic_gene2GOs:
		dic_gene2GOs[gene] = set()
	dic_gene2GOs[gene].add(GO)

if args.outfmt == 'gene2GO' or args.outfmt == 'gene2GOname':
	for gene, set_GO in sorted(dic_gene2GOs.items()):
		for GO in set_GO:
			if args.outfmt == 'gene2GO':
				OF.write(gene+'\t'+GO+'\n')
			elif args.outfmt == 'gene2GOname':
				GOname = GO+'+'+dic_GO2name[GO]
				OF.write(gene+'\t'+GOname+'\n')
elif args.outfmt == 'gene2GOs' or args.outfmt == 'gene2GOnames':
	for gene, set_GO in sorted(dic_gene2GOs.items()):
		if args.outfmt == 'gene2GOs':
			OF.write(gene+'\t'+','.join(sorted(set_GO))+'\n')
		elif args.outfmt == 'gene2GOnames':
			lst_GOnames=map(lambda x: x+'+'+dic_GO2name[x], sorted(set_GO))
			OF.write(gene+'\t'+','.join(lst_GOnames)+'\n')

elif args.outfmt == 'GO2genes' or args.outfmt == 'GOname2genes':
	for GO, set_gene in sorted(dic_GO2genes.items()):
		if args.outfmt == 'GO2genes':
			OF.write(GO+'\t'+','.join(sorted(set_gene))+'\n')
		elif args.outfmt == 'GOname2genes':
			GOname = GO+'+'+dic_GO2name[GO]
			OF.write(GOname+'\t'+','.join(sorted(set_gene))+'\n')
