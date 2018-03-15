#!/usr/bin/env python3

# JA script for [GO|KOG|.etc.] term category enrichment over a set of clusters

from scipy.stats import hypergeom as hg
import pandas as pd
import sys,math
import pickle as p

def enrichment(ff,grps,idkey,termkey,termtype,termname,typefilter):
	t = pd.read_table(ff,sep='\t',dtype='str')
	result = {}
	termnames = {}
	pcut = 0.05
	for i,row in t.iterrows():
		if not row[termkey] in termnames: termnames[row[termkey]] = row[termname]
	# for filtering lists with multiple combined categories (e.g. old JGI GO table)
	if not typefilter=='':
		t = t[ t[typefilter].map(lambda x: x==termtype) ]
	ll = list(t[termkey])
	terms = set(ll)
	# bg counts
	bg = {}
	for term in terms:
		bg[term] = ll.count(term)
	# for each group of ids
	for k,g in grps:
		# subtable for only this group of ids
		tg = t[ t[idkey].map(lambda x: x in g) ]
		# this table might be empty if no terms
		if tg.shape[0]==0: continue
		llg = list(tg[termkey])
		#, for each term, compare group counts to background counts
		gterms = []
		for term in terms:
			nt_g = llg.count(term)
			p = hg(t.shape[0],bg[term],tg.shape[0]).pmf(nt_g)
			if p>pcut: continue
#			print(t.shape[0], bg[term], tg.shape[0], nt_g, p)
#			return
			gterms.append((term,termnames[term],p))
		if len(gterms)>0: result[k] = gterms
	return result

info = sys.argv[1]
clsf = sys.argv[2]
key = 'cl'
if len(sys.argv)>3: key = sys.argv[3]

# for Thaps, Nanoc, Fracy, Psemu, Phatr descriptive GO terms (old JGI)
idkey = 'proteinId'
termkey = 'goAcc'
typefilter = 'gotermType'
termname = 'goName'
termtype = 'biological_process'

## for Cyccr KOG
## this needs to be updated with InterProScan 
#idkey = 'augustus_id'
#termkey = 'react'
#typefilter = ''
#termname = 'react'
#termtype = ''
#
## for Fracy, Nanoc, Psemu KOG
#idkey = 'proteinId'
#termkey = 'kogClass'
#typefilter = ''
#termname = 'kogClass'
#termtype = ''

cls = {}
for i,row in pd.read_table(clsf,dtype='str').iterrows():
	clid = str(row[key])
	if clid=='NA' or clid=='nan': continue
	if not clid in cls: cls[clid] = []
	cls[clid].append(row['id'])
print('%i cls loaded with %i unique keys'%(len(cls),len(set(cls.keys()))))

e = enrichment(info,cls.items(),idkey,termkey,termtype,termname,typefilter)
p.dump(e,open('terms.%s.p'%key,'wb'))
