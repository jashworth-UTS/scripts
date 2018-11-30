#!/usr/bin/env python
# remedial script for calculating restriction enzyme digest band sizes for CIRCULAR (plasmid) sequences
# more/different code would be needed for LINEAR fragment digests
# also computes band sizes corresponding to intermediate fragments containing one uncut site

import re, os, sys
from itertools import combinations as cmb

class Site:
	def __init__(self,name,seq):
		self.name = name
		self.seq = seq
	def __str__(self):
		return('%s: %s' %(self.name, self.seq))

sites = [
	Site('NheI','GCTAGC'),
	Site('EcoRI','GAATTC'),
]

# a 'regular expression' object for matching the restriction site occurrences
if len(sites) == 1:
	re_anysite = re.compile(sites[0].seq,re.IGNORECASE)
else:
	re_anysite = re.compile('|'.join([s.seq for s in sites]),re.IGNORECASE)

def readfasta(fname):
	totread = 0
	seqs = []
	seqname=''
	seq=[]
	sys.stderr.write('reading sequence from %s\n' %fname)
	for line in open(fname):
		if line.startswith('>'):
			# fasta name line
			if seqname != '' and len(seq)>0:
				# if we already haves some sequence loaded, assume end of previous seq, store
				seqs.append( (seqname,''.join(seq)) )
				totread += len(seqs[-1][1])
			# store new seqname in name line
			seqname=re.sub('>','',line.strip()).strip()
			# start new sequence
			seq=[]
		else:
			# is not a name line, so add line to current sequence lines
			seq.append( re.sub(' ','', line.strip() ) )
	# last hanging seq
	if seqname != '' and len(seq)>0:
		seqs.append( (seqname,''.join(seq)) )
		totread += len(seqs[-1][1])
	sys.stderr.write('%i letters read\n' %totread)
	return(seqs)

searchseqs = readfasta(sys.argv[1])

# recursive function that adds up every number of possible fragment combinations in a nested way
def do_combos(elems,output,n):
	if n==0: return
	for celems in cmb(elems,n):
		output.append( sum(celems) )
	do_combos(elems,output,n-1)

for searchseq in searchseqs:
	print('Sequence: %s (length: %i)' %(searchseq[0], len(searchseq[1])))

	# find site locations
	locs = []
	for m in re_anysite.finditer(searchseq[1]): locs.append(m.start())
	sys.stderr.write('site locations: %s\n' %str(locs))

	# adjust site locations (set first occurrence to zero)
	locs2 = []
	for l in locs: locs2.append(l-locs[0])
	sys.stderr.write('reindexed locations: %s\n' %str(locs2))
	# add full linear length of sequence for band size calcs
	locs2.append(len(searchseq[1]))

	# compute band sizes
	sizes = []
	for i in range(len(locs2)-1):
		sizes.append( locs2[i+1] - locs2[i] )
	sys.stderr.write('sizes of all-cut fragments: %s\n' %str(sizes))

	# compute all possible band sizes including incomplete digestion
	allsizes = []
	do_combos(sizes,allsizes,len(sizes))

	sys.stderr.write('sizes of all calculated fragments: %s\n' %str(allsizes))

	allsizes.sort()
	allsizes = set(allsizes)
	print('Sizes of all fragments calculated:')
	for s in sorted(allsizes): print('\t%i' %s)
