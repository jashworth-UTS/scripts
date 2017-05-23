#!/usr/bin/env python

import sqlite3 as sq
con = sq.connect('diatoms.sqlite.db')
c = con.cursor()

from reciprocal_dict import Model

import pickle
def loadhoms(file='reciprocal_dict.p'):
	return pickle.load(open(file,'rb'))
# do this manually--slow
#hd = loadhoms()

def write_fasta(seqs,fname):
	f=open(fname,'w')
	for name,seq in seqs:
		f.write('>%s\n%s\n' %(name,seq))
	f.close()

def thapsgrep(pattern):
#	import pandas as pd
	if not ann: ann = pd.read_csv('/Users/justin/diatoms/genemodels/tps.all.models.ja.new.tsv',delimiter='\t')
	cnd = ann['name'].str.contains(pattern, na=False)
	ids = list( ann['id'][cnd] )
	return [ Model('Thaps',id) for id in ids ]

def orgid(shortname):
	return c.execute('select id from org where short="%s"' %shortname).fetchone()[0]

def cdsid(model):
	return c.execute('select id from cds where org=%s AND pid="%s"' %(orgid(model.org),model.id)).fetchone()[0]

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

# ...map orthologous/homologous genes from reciprocal dicts into db...
# ...smart/greedy extension of putative real CDS's and resulting plausible true promoter coords by protein alignment across homologs
# ...write promoters for sets of sets of orthologous/homologous genes...
# Lhcs
# ACT1's
# EFs
# metabolic pathways
# carbon concentrating?
# anything with a cool signal in the expression clustering
