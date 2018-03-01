#!/usr/bin/env python

import sys,re
from functools import total_ordering

''' JGI GFF example
scaffold_1	JGI	exon	2702	2844	.	-	.	name "fgenesh_newKGs_kg.1__1__EST_ALL.fasta.Contig5314"; transcriptId 349042
scaffold_1	JGI	exon	2995	3157	.	-	.	name "fgenesh_newKGs_kg.1__1__EST_ALL.fasta.Contig5314"; transcriptId 349042
scaffold_1	JGI	CDS	3042	3157	.	-	0	name "fgenesh_newKGs_kg.1__1__EST_ALL.fasta.Contig5314"; proteinId 349042; exonNumber 3
scaffold_1	JGI	stop_codon	3042	3044	.	-	0	name "fgenesh_newKGs_kg.1__1__EST_ALL.fasta.Contig5314"
scaffold_1	JGI	exon	3042	3176	.	-	.	name "fgenesh_newKGs_kg.1__1__EST_ALL.fasta.Contig5314"; transcriptId 349042
scaffold_1	JGI	CDS	3042	3176	.	-	2	name "fgenesh_newKGs_kg.1__1__EST_ALL.fasta.Contig5314"; proteinId 349042; exonNumber 2
scaffold_1	JGI	stop_codon	3042	3044	.	-	0	name "fgenesh_newKGs_kg.1__1__EST_ALL.fasta.Contig5314"
scaffold_1	JGI	exon	3442	3678	.	-	.	name "fgenesh_newKGs_kg.1__1__EST_ALL.fasta.Contig5314"; transcriptId 349042
scaffold_1	JGI	CDS	3442	3502	.	-	2	name "fgenesh_newKGs_kg.1__1__EST_ALL.fasta.Contig5314"; proteinId 349042; exonNumber 1
scaffold_1	JGI	start_codon	3500	3502	.	-	0	name "fgenesh_newKGs_kg.1__1__EST_ALL.fasta.Contig5314"
'''

''' what Stringtie wants (GFF3)
ctg123	example	mRNA	1300	9950	.	+	.	ID=t_012143;gene_name=EDEN
ctg123	example	exon	1300	1500	.	+	.	Parent=t_012143
ctg123	example	exon	3000	3902	.	+	.	Parent=t_012143
ctg123	example	CDS	3301	3902	.	+	0	Parent=t_012143
ctg123	example	exon	5000	5500	.	+	.	Parent=t_012143
ctg123	example	CDS	5000	5500	.	+	1	Parent=t_012143
ctg123	example	exon	7000	9000	.	+	.	Parent=t_012143
ctg123	example	CDS	7000	7600	.	+	1	Parent=t_012143
ctg123	example	exon	9400	9950	.	+	.	Parent=t_012143
'''

skipkeys = [
	'name'
]

@total_ordering
class Feature:
	sep = '\t'
	def __init__(self,seq=None,source=None,type=None,start=None,end=None,score='.',strand='.',phase='.',attributes={}):
		self.seq    = seq
		self.source = source
		self.type   = type
		self.start  = start
		self.end    = end
		self.score  = score
		self.strand = strand
		self.phase  = phase
		self.attributes = attributes
	def __str__(self):
		attstr = ';'.join(['%s=%s' %(k,v) for k,v in self.attributes.items()])
		return self.sep.join( [self.seq, self.source, self.type, str(self.start), str(self.end), self.score, self.strand, self.phase, attstr])
	def __eq__(self,other):
		return self == other
	def __ne__(self,other):
		return self != other
	def __lt__(self,other):
		if self.seq < other.seq: return True
		if self.start < other.start: return True
		if self.end < other.end: return True
		return False

class Parent(Feature):
	def __init__(self,id):
		Feature.__init__(self)
		self.id = id
		self.features = []
		self.attributes = {'ID':self.id}
	def summarize(self):
		self.start = -1
		self.end = -1
		self.strand = None
		self.seq = None
		self.source = None
		self.type = 'mRNA'
		for f in self.features:
			if self.strand==None: self.strand = f.strand
			elif self.strand != f.strand: sys.stderr.write('WARNING: child feature strand mismatch!! %s\n' %str(f))
			if self.seq==None: self.seq = f.seq
			if self.source==None: self.source = f.source
			if self.start == -1: self.start = f.start
			if self.end == -1: self.end = f.end
			if f.start < self.start: self.start = f.start
			if f.end > self.end: self.end = f.end
			f.attributes['Parent'] = self.id
			for pk in parent_keys:
				if pk in f.attributes: f.attributes.pop(pk)

parents_by_contig = {}

parent_keys = [
	'transcriptId',
	'transcript_id',
]

for line in sys.stdin:
	l = line.strip().split('\t')
	atts = {}
	for att in l[8].split(';'):
		att = re.split('[ =]',att.strip())
		k = att[0]
		if k in skipkeys: continue
		atts[k] = att[1]
	contig = l[0]
	if not contig in parents_by_contig: parents_by_contig[contig] = {}
	f = Feature(l[0], l[1], l[2], int(l[3]), int(l[4]), l[5], l[6], l[7], atts)

	for pk in parent_keys:
		if pk in atts:
			id = atts[pk]
			if not id in parents_by_contig[contig]: parents_by_contig[contig][id] = Parent(id)
			parents_by_contig[contig][id].features.append(f)

for c in sorted(parents_by_contig):
	for k in parents_by_contig[c]:
		parents_by_contig[c][k].summarize()

for c in sorted(parents_by_contig):
	for p in sorted(parents_by_contig[c].values()):
		print(str(p))
		for f in p.features:
			print(str(f))
