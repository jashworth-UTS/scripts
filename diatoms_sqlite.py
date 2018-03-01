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

def annotation(model):
	ann = (domain(model), superfam(model))
	return(ann)

def domain(model):
	oid = get_orgid(model.org)
#	print('fetching domain for org %s (id %i) pid %s' %(model.org, oid, model.id))
	return c.execute('select * from domain_hits where org=%i and id="%s"' %(oid,model.id)).fetchall()

def superfam(model):
	oid = get_orgid(model.org)
	return c.execute('select * from superfam_hits where org=%i and id="%s"' %(oid,model.id)).fetchall()

def writefasta(seqs,fname):
	f=open(fname,'w')
	for s in seqs:
		if not s: continue
		f.write('>%s\n%s\n' %(s[0],s[1]))
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
		return c.execute('select id from cds where org=%s AND pid="%s"' %(get_orgid(model.org),model.id)).fetchone()[0]
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

#def update_noncoding_boundaries():
	# systematically determine and update coordinates for 5' and 3' "noncoding regions" adjacent to models (as defined by adjacent model boundaries) in DB table

def modelseq(orgid,seq,start,end,strand):
	s = c.execute('select seq from genomic where org = %i AND name = "%s"' %(orgid,seq)).fetchone()[0]
	l = len(s)
	# note that genome is 1-indexed, and Python is 0-indexed, and Python list/str indexing is [,)
	s = s[(start-1):end]

	# 1 is antisense
	if strand == 1: return rvs_comp_str(s)
	return s

# upstream seq relative to one cds or gene start (by internal db id)
def upstreamseq(orgid,seq,start,end,strand,us=500,ds=100):
	s = c.execute('select seq from genomic where org = %i AND name = "%s"' %(orgid,seq)).fetchone()[0]
	l = len(s)
	# note that genome is 1-indexed, and Python is 0-indexed
	start = start-1
	end = end-1

	# 1 is antisense
	if strand == 1:
		left = end - ds
		right = end + us
		# bounds correction
		if left < 1: left = 1
		if right > l: right = l
		# note that 'end' is last included position, requiring a register shift of Python list/str indexing
		s = s[(left+1):(end+1)].upper() + s[(end+1):(right+1)].lower()
		return rvs_comp_str(s)

	left = start - us
	right = start + ds
	# bounds correction
	if left < 1: left = 1
	if right > l: right = l
	return s[left:start].lower() + s[start:right].upper()

# downstream seq relative to a cds or gene end (by internal db id)
def downstreamseq(orgid,seq,start,end,strand,us=0,ds=500):
	s = c.execute('select seq from genomic where org = %i AND name = "%s"' %(orgid,seq)).fetchone()[0]
	l = len(s)
	# note that genome is 1-indexed, and Python is 0-indexed
	start = start-1
	end = end-1

	# 1 is antisense
	if strand == 1:
		left = start - ds
		right = start + us
		# bounds correction
		if left < 1: left = 1
		if right > l: right = l
		s = s[left:start].lower() + s[start:right].upper()
		return rvs_comp_str(s)

	left = end - us
	right = end + ds
	# bounds correction
	if left < 1: left = 1
	if right > l: right = l
	# note that 'end' is last included position, requiring register shift of Python list/str indexing
	return s[(left+1):(end+1)].upper() + s[(end+1):(right+1)].lower()

def get_model_coords(model,adjacent_noncoding_boundaries=False):
	coords = None
	orgid = get_orgid(model.org)
	# see if model id is in cds
	result = c.execute('select seq,start,end,strand from cds where org=%s and pid="%s"' %(orgid,model.id))
	try:
		(seq,start,end,strand) = result.fetchone()
		if adjacent_noncoding_boundaries:
			prev_end = c.execute('select max(end) from cds where org=%s and seq="%s" and end < %s' %(orgid,seq,start)).fetchone()[0]
			next_start = c.execute('select min(start) from cds where org=%s and seq="%s" and start > %s' %(orgid,seq,end)).fetchone()[0]
			coords = (orgid,seq,start,end,strand,prev_end,next_start)
		else: coords = (orgid,seq,start,end,strand)
		print('cds:',model.org,model.id,coords)
	except:
		# see if model id is in gene
		result = c.execute('select seq,start,end,strand from gene where org=%s and gid="%s"' %(orgid,model.id))
		try:
			(seq,start,end,strand) = result.fetchone()
			if adjacent_noncoding_boundaries:
				prev_end = c.execute('select max(end) from gene where org=%s and seq="%s" and end < %s' %(orgid,seq,start)).fetchone()[0]
				next_start = c.execute('select min(start) from gene where org=%s and seq="%s" and start > %s' %(orgid,seq,end)).fetchone()[0]
				coords = (orgid,seq,start,end,strand,prev_end,next_start)
			else: coords = (orgid,seq,start,end,strand)
			print('gene:',model.org,model.id,coords)
		except:
			print('lookup failed for org %s (id %i) model id %s' %(model.org,orgid,model.id))
			return None
	return coords

def seq_for_model(model):
	coords = get_model_coords(model)
	if coords: (orgid,seq,start,end,strand) = coords
	else: return coords
	return ('%s_%s' %(model.org,model.id), modelseq(orgid,seq,start,end,strand))

def promoter_for_model(model,trim=True,us=500,ds=0):
	coords = get_model_coords(model,adjacent_noncoding_boundaries=True)
	if coords: (orgid,seq,start,end,strand,prev_end,next_start) = coords
	else: return coords
	# trim off preceding gene model sequence (return non-coding only)
	if trim:
		if strand == 1:
			if next_start: us = min(us,next_start-1-end)
		else:
			if prev_end: us = min(us,start-prev_end+1)
	return ('%s_%s' %(model.org,model.id), upstreamseq(orgid,seq,start,end,strand,us=us,ds=ds))

def threeprime_for_model(model,trim=True,us=0,ds=500):
	coords = get_model_coords(model,adjacent_noncoding_boundaries=True)
	if coords: (orgid,seq,start,end,strand,prev_end,next_start) = coords
	else: return coords
	# trim off following gene model sequence (return non-coding only)
	if trim:
		if strand == 1:
			if prev_end: ds = min(ds,start-prev_end+1)
		else:
			if next_start: ds = min(ds,next_start-1-end)
	return ('%s_%s' %(model.org,model.id), downstreamseq(orgid,seq,start,end,strand,us=us,ds=ds))

# ...map orthologous/homologous genes from reciprocal dicts into db...
# ...smart/greedy extension of putative real CDS's and resulting plausible true promoter coords by protein alignment across homologs
# ...write promoters for sets of sets of orthologous/homologous genes...
# Lhcs
# ACT1's
# EFs
# metabolic pathways
# carbon concentrating?
# anything with a cool signal in the expression clustering
