#!/usr/bin/env python

import sqlite3 as sq
con = sq.connect('diatoms.sqlite.db')
c = con.cursor()

from reciprocal_dict import *
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

def seq_from_coords(orgid,seq,start,end,strand):
	s = c.execute('select seq from genomic where org = %i AND name = "%s"' %(orgid,seq)).fetchone()[0]
	l = len(s)
	# note that genome is 1-indexed, and Python is 0-indexed
	s = s[(start-1):end]

	# 1 is antisense
	if strand == 1: return rvs_comp_str(s)
	return s

# promoter seq region relative to one cds (by cds internal db id)
def promseq(orgid,seq,start,end,strand,us=500,ds=100):
	s = c.execute('select seq from genomic where org = %i AND name = "%s"' %(orgid,seq)).fetchone()[0]
	l = len(s)
	# note that genome is 1-indexed, and Python is 0-indexed
	start = start-1
	end = end-1

	# 1 is antisense
	if strand == 1:
		left = end - ds+1
		right = end + us
		# bounds correction
		if left < 1: left = 1
		if right > l: right = l
		# note that 'end' is last included position
		s = s[left:(end+1)] + s[(end+1):right].lower()
		return rvs_comp_str(s)

	left = start - us+1
	right = start + ds
	# bounds correction
	if left < 1: left = 1
	if right > l: right = l
	return s[left:start].lower() + s[start:right]

def get_model_coords(model):
	orgid = get_orgid(model.org)
	# see if model id is in cds
	result = c.execute('select seq,start,end,strand from cds where org=%s and pid="%s"' %(orgid,model.id))
	try:
		(seq,start,end,strand) = result.fetchone()
		print('cds:',model.org,orgid,seq,start,end,strand)
	except:
		# see if model id is in gene
		result = c.execute('select seq,start,end,strand from gene where org=%s and gid="%s"' %(orgid,model.id))
		try:
			(seq,start,end,strand) = result.fetchone()
			print('gene:',model.org,orgid,seq,start,end,strand)
		except:
			print('lookup failed for org %s (id %i) model id %s' %(model.org,orgid,model.id))
			return None
	return (orgid,seq,start,end,strand)

def seq_for_model(model):
	coords = get_model_coords(model)
	if coords: (orgid,seq,start,end,strand) = coords
	else: return coords
	return ('%s_%s' %(model.org,model.id), seq_from_coords(orgid,seq,start,end,strand))

def promoter_for_model(model,**kwargs):
	coords = get_model_coords(model)
	if coords: (orgid,seq,start,end,strand) = coords
	else: return coords
	return ('%s_%s' %(model.org,model.id), promseq(orgid,seq,start,end,strand,**kwargs))

# ...map orthologous/homologous genes from reciprocal dicts into db...
# ...smart/greedy extension of putative real CDS's and resulting plausible true promoter coords by protein alignment across homologs
# ...write promoters for sets of sets of orthologous/homologous genes...
# Lhcs
# ACT1's
# EFs
# metabolic pathways
# carbon concentrating?
# anything with a cool signal in the expression clustering
