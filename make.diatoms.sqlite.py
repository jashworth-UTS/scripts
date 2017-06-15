#!/usr/bin/env python
import os,re
import Fasta

# PURPOSE: build a convenient and orderly sqlite3 database for multi-species record-keeping, lookups, references to genome features, etc for multi-species comparative work
# This script builds the database from [properly formatted and named] genome files
# Then interact interactive or progammatically with database in python or ipython via sqlite3 API
# see diatoms.sqlite.py for usage functions

import sqlite3 as sq
con = sq.connect('diatoms.sqlite.db')
c = con.cursor()

# files needed for this to build a DB:
# [org].gff (organism genome feature models in gff/gff3 format) /at least/ CDS features that have 'protein_id=##' specified in the attributes (ninth) column
# [org].aa.ids.fa (protein sequences in fasta files with *only* protein_id as fasta header

orgs = [
	"Arath",
	"Chlre",
	"Cyccr",
	"Emihu",
	"Fracy",
	"Phatr2",
	"Phatr3",
	"Physo",
	"Psemu",
	"Thaoc",
	"Thaps",
]

dropall = True
#dropall = False
if dropall:
	print('dropping old tables')
	# try/pass is silent if old tables do not exist
	try: c.execute('drop table org')
	except: pass
	try: c.execute('drop table genomic')
	except: pass
	try: c.execute('drop table cds')
	except: pass
	try: c.execute('drop table gene')
	except: pass
	try: c.execute('drop table pep')
	except: pass
	con.commit()

print('creating tables')
tables = [t[0] for t in c.execute('SELECT name FROM sqlite_master WHERE type="table"').fetchall()]
if not 'org' in tables: c.execute('create table org (id integer primary key autoincrement, short text)')
if not 'genomic' in tables: c.execute('create table genomic (id integer primary key autoincrement, org integer, name text, seq text)')
if not 'cds' in tables: c.execute('create table cds (id integer primary key autoincrement, org integer, seq text, pid text, start integer, end integer, strand boolean)')
if not 'gene' in tables: c.execute('create table gene (id integer primary key autoincrement, org integer, seq text, gid text, start integer, end integer, strand boolean)')
if not 'pep' in tables: c.execute('create table pep (id integer primary key autoincrement, org integer, name text, seq text)')
con.commit()

#print( c.execute('select * from org').fetchall() )

if not dropall:
	# didn't drop tables for complete replacement, delete/replace by org
	# (faster for adding/refreshing data for just some orgs)
	for org in orgs:
		orgid = c.execute('select id from org where short="%s"' %org).fetchall()
		print(orgid)
		con.commit()
		for id in orgid:
			print('removing org %s (%i) from tables' %(org,id[0]))
			c.execute('delete from org where id=%s' %id[0])
			c.execute('delete from genomic where org=%s' %id[0])
			c.execute('delete from cds where org=%s' %id[0])
			c.execute('delete from gene where org=%s' %id[0])
			c.execute('delete from pep where org=%s' %id[0])
			con.commit()

class Feature:
	def __init__(self,org,seq,fid,start,end,strand):
		self.org    = org
		self.seq    = seq
		self.fid    = fid
		self.start  = int(start)
		self.end    = int(end)
		self.strand = strand

	def extend(self,coords):
		for i in coords:
			if self.start > int(i): self.start = int(i)
			if self.end < int(i): self.end = int(i)

def read_gff(gff):

	features = {
		'CDS' : {},
		'gene' : {},
	}

	if not os.path.exists(gff):
		print('file %s not found' %gff)
		return features

	for l in open(gff):
		if l.startswith('#'): continue
		s = l.strip().split('\t')
		featuretype = s[2]
		if not featuretype in features: continue

		strand = 0
		if s[6] == '-': strand = 1

		featureid = ''
		atts = s[8].replace(',',';').split(';')
		for att in atts:
			att = att.strip()
			kv = att.split('=')
			if featuretype == 'CDS' and kv[0] == 'protein_id': featureid = kv[1]
			elif featuretype == 'gene' and kv[0] == 'gene_id': featureid = kv[1]
		if featureid == '': continue

		start = int(s[3])
		end = int(s[4])
		if not featureid in features[featuretype]:
			features[featuretype][featureid] = Feature(orgid,s[0],featureid,start,end,strand)
		else:
			features[featuretype][featureid].extend([start,end])

	for ft in features:
		print('read %i %s features from %s' %(len(features[ft]),ft,gff))
	return features

### MAIN ###
### fill tables by org
for org in orgs:
	c.execute('insert into org (short) values ("%s")' %org)
	con.commit()
	orgid = c.execute('select id from org where short="%s"' %org).fetchone()[0]

	print('%i %s' %(orgid,org))

	# genomic contig sequences
	gen = '%s.genome.fa' %org
	if not os.path.exists(gen): continue
	gs = Fasta.FastaSeqs()
	gs.loadseqs([gen])
	print('%i genomic contig seqs' %len(gs.seqs))
	for s in gs.seqs.values():
		ex = 'insert into genomic(org,name,seq) values(%i,"%s","%s")' %(orgid,s.name,s.seq)
		c.execute(ex)
	con.commit()

	gff = '%s.gff' %org
	features = read_gff(gff)

	# CDS models (seems to be present and usable in all gffs)
	for cds in features['CDS'].values():
		ex = 'insert into cds (org,seq,pid,start,end,strand) values (%i,"%s","%s","%s","%s","%s")' %(orgid,cds.seq,cds.fid,cds.start,cds.end,cds.strand)
		c.execute(ex)
	con.commit()

	# gene models (missing currently in majority of diatom gffs....)
	for gene in features['gene'].values():
		ex = 'insert into gene (org,seq,gid,start,end,strand) values (%i,"%s","%s","%s","%s","%s")' %(orgid,gene.seq,gene.fid,gene.start,gene.end,gene.strand)
		c.execute(ex)
	con.commit()

	# peptide sequences
	pep = '%s.aa.ids.fa' %org
	if not os.path.exists(pep): continue
	ss = Fasta.FastaSeqs()
	ss.loadseqs([pep])
	print('%i pep seqs' %len(ss.seqs))
	for s in ss.seqs.values():
		ex = 'insert into pep(org,name,seq) values(%i,"%s","%s")' %(orgid,s.name,s.seq)
		c.execute(ex)
	con.commit()

c.execute('select id,short from org')
orgs = c.fetchall()
print('contents of orgs table:')
print(orgs)

con.close()

# see diatoms.sqlite.py for usage functions
