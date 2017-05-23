#!/usr/bin/env python

import sqlite3 as sq
con = sq.connect('diatoms.sqlite.db')
c = con.cursor()

from reciprocal_dict import Model
from AshworthUtil import *

import pickle
def loadhoms(file='reciprocal_dict.p'):
	return pickle.load(open(file,'rb'))
# do this manually--slow
#hd = loadhoms()

def writefasta(seqs,fname):
	f=open(fname,'w')
	for name,seq in seqs:
		f.write('>%s\n%s\n' %(name,seq))
	f.close()

def thapsgrep(pattern):
	import pandas as pd
	ann = pd.read_csv('/Users/justin/diatoms/genemodels/tps.all.models.ja.new.tsv',delimiter='\t')
	cnd = ann['name'].str.contains(pattern, na=False)
	ids = list( ann['id'][cnd] )
	return [ Model('Thaps',id) for id in ids ]

def get_orgid(shortname):
	return c.execute('select id from org where short="%s"' %shortname).fetchone()[0]

def get_cdsid(model):
	try:
		return c.execute('select id from cds where org=%s AND pid="%s"' %(get_orgid(model.org),model.id)).fetchone()
	except:
		pass

def bad_coordinates(table):
	# check all(starts<=ends in loaded gff models)
	starts = c.execute('select start from %s' %table).fetchall()
	ends = c.execute('select end from %s' %table).fetchall()
	con.commit()
	lt = [a<=b for (a,b) in zip(starts,ends)]
	f = [i for i in range(len(r)) if not r[i]]
	for i in f:
		print( c.execute('select * from %s where id=%i'%(table,i+1)).fetchall())
	return(f)

# promoter seq region relative to one cds (by cds internal db id)
def promseq(cds,upstream=500,downstream=100):
	if not cds: return ''
	(id,org,seq,pid,start,end,strand) = c.execute('select * from cds where id = %i' %cds).fetchone()
	ps = start - upstream
	pe = end + downstream
	if strand == 1:
		ps = start - downstream
		pe = end + upstream

	if ps < 1: ps = 1

	s = c.execute('select seq from genomic where name = "%s" AND org = "%s"' %(seq,org)).fetchone()[0]
	l = len(s)
	if ps >= l:
		print('out of bounds, check feature %s and sequence %s for org %s' %(cds,seq,org))
	if pe > l: pe = l
	s = s[ps-1:pe]
	if strand == 1:
		s = rvs_comp_str(s)
	return(s)

def promoter_for_model(model):
	return ('%s_%s' %(model.org,model.id), promseq(get_cdsid(model)))

# ...map orthologous/homologous genes from reciprocal dicts into db...
# ...smart/greedy extension of putative real CDS's and resulting plausible true promoter coords by protein alignment across homologs
# ...write promoters for sets of sets of orthologous/homologous genes...
# Lhcs
# ACT1's
# EFs
# metabolic pathways
# carbon concentrating?
# anything with a cool signal in the expression clustering
