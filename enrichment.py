#!/usr/bin/env python3

# JA script for [GO|KOG|.etc.] term category enrichment over a set of clusters

from scipy.stats import hypergeom as hg
import pandas as pd
import sys,math
import pickle as p

def enrichment(infofile,grps,termidkey,termkey,typefilter,termtype,termname):
	t = pd.read_table(infofile,dtype='str')
	result = {}
	termnames = {} # long description for terms
	pcut = 0.01
	for i,row in t.iterrows():
		if not row[termkey] in termnames: termnames[row[termkey]] = row[termname]
	# for filtering lists with multiple combined categories (e.g. old JGI GO table)
	if not typefilter=='':
		t = t[ t[typefilter].map(lambda x: x==termtype) ]
	allterms = list(t[termkey]) # full column of term keys (e.g. term accession numbers)
	allids = list(t[termidkey]) # full column of protein ids
	# bg counts
	bg = {}
	uniqueterms = set(allterms) # unique set of all term keys
	for term in uniqueterms: bg[term] = allterms.count(term) # background counts for each term
	# for each group of ids
	for clkey,grp in grps:
		# subtable for only this group of ids
		rows = [] # table/list indices corresponding to the ids in this group
		for i in range(len(allids)):
			if allids[i] in grp: rows.append(i)
		if len(rows)==0: continue
		groupterms = [allterms[r] for r in rows] # column of term keys for this group of ids
		#, for each term, compare group counts to background counts
		gterms = []
		for term in set(groupterms): # for all terms observed for this group
			p = hg(len(allterms),bg[term],len(groupterms)).pmf(groupterms.count(term))
			if p>pcut: continue
			gterms.append((term,termnames[term],p))
		if len(gterms)>0: result[clkey] = gterms
	return result

infofile = sys.argv[1]
clsf = sys.argv[2]
clkey = 'cl'
if len(sys.argv)>3: clkey = sys.argv[3]
clidkey = 'id'
if len(sys.argv)>4: clidkey = sys.argv[4]

# for Thaps, Phatr, Emihu, Nanoc, Fracy, Psemu, Phatr descriptive GO terms (JGI style)
termidkey = 'proteinId'
termkey = 'goAcc'
typefilter = ''
# uncomment and set termtype to limit to subtype of GO terms
#typefilter = 'gotermType'
termname = 'goName'
termtype = 'biological_process'

## for Cyccr KOG
# needs to be updated with InterProScan to get goterms
if 'reactome' in infofile:
	termidkey = 'proteinId'
	termkey = 'accession'
	typefilter = ''
	termname = 'reactome'
	termtype = ''

## for Fracy, Nanoc, Psemu KOG
#termidkey = 'proteinId'
#termkey = 'kogClass'
#typefilter = ''
#termname = 'kogClass'
#termtype = ''

cls = {}
for i,row in pd.read_table(clsf,dtype='str').iterrows():
	clid = str(row[clkey])
	if clid=='NA' or clid=='nan': continue
	if not clid in cls: cls[clid] = []
	cls[clid].append(row[clidkey])
# ensure non-redundant cluster member ids
for c,vs in cls.items(): cls[c] = set(vs)

print('%i groups loaded, infofile %s'%(len(cls),infofile))

e = enrichment(infofile,cls.items(),termidkey,termkey,typefilter,termtype,termname)
p.dump(e,open('terms.%s.p'%clkey,'wb'))
