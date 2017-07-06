#!/usr/bin/env python
import sys

# PURPOSE: add conserved domain database (CDD) hits from rpsblast/rpsbproc to diatom database

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

cdd_tables = ['domain_hits','superfam_hits','motif_hits','site_hits']

# set up tables

dropall = True
#dropall = False
if dropall:
	print('dropping old tables')
	for t in cdd_tables:
		# try/pass is silent if old tables do not exist
		try:
			c.execute('drop table %s' %t)
			print('dropping old table %s' %t)
		except: pass
	con.commit()

tables = [t[0] for t in c.execute('SELECT name FROM sqlite_master WHERE type="table"').fetchall()]

if not 'domain_hits' in tables:
	print('creating table domain_hits')
	c.execute('create table domain_hits (id integer primary key autoincrement, org integer, q text, hittype text, pssmid integer, hit_from integer, hit_to integer, eval real, bitscore real, accession text, shortname text, desc text)')
con.commit()

if not 'superfam_hits' in tables:
	print('creating table superfam_hits')
	c.execute('create table superfam_hits (id integer primary key autoincrement, org integer, q text, hittype text, pssmid integer, hit_from integer, hit_to integer, eval real, bitscore real, accession text, shortname text, desc text)')
con.commit()

#print( c.execute('select * from org').fetchall() )

files = sys.argv[1:]

eval_cut = 1e-6

for f in files:
	org = f.strip().split('.')[0]
	if not org in orgs:
		print('file/org %s/%s not in orgs list--check file against orgs list' %(f,org))
		continue

	# first delete old records in CDD table(s) for this org
	orgid = c.execute('select id from org where short="%s"' %org).fetchone()[0]
	for t in cdd_tables:
		if not t in tables: continue
		print('removing org %s (%i) from %s table' %(org,orgid,t))
		c.execute('delete from %s where org=%s' %(t,orgid))
		con.commit()

	# parse CDD rpsblast/rpsbproc file, add records
	domain_hits = []
	superfam_hits = []
	print('parsing input...')
	for l in open(f):
		if l.startswith('DOMAINS'):
			s = l.strip().split('\t')
			ev = float(s[8])
			if ev<=eval_cut:
				domain_hits.append( (orgid,s[3],s[4],s[5],s[6],s[7],s[8],s[9],s[10],s[11],s[14]) )
		elif l.startswith('SUPERFAMILIES'):
			s = l.strip().split('\t')
			ev = float(s[8])
			if ev<=eval_cut:
				superfam_hits.append( (orgid,s[3],s[4],s[5],s[6],s[7],s[8],s[9],s[10],s[11],s[14]) )

	print('inserting into DB...')
	c.executemany('insert into domain_hits(org,q,hittype,pssmid,hit_from,hit_to,eval,bitscore,accession,shortname,desc) values(?,?,?,?,?,?,?,?,?,?,?)', domain_hits)
	con.commit()
	print('added %i domain hits to cdd domain_hits table for org %s (id %i)' %(len(domain_hits),org,orgid))

	c.executemany('insert into superfam_hits(org,q,hittype,pssmid,hit_from,hit_to,eval,bitscore,accession,shortname,desc) values(?,?,?,?,?,?,?,?,?,?,?)', superfam_hits)
	con.commit()
	print('added %i superfamily hits to cdd superfam_hits table for org %s (id %i)' %(len(superfam_hits),org,orgid))

con.close()

# see diatoms.sqlite.py for usage functions
