#!/usr/bin/env python

import Fasta
from AshworthUtil import *
import sys,re

primer = sys.argv[1]
f = sys.argv[2]

def expandprimer(primer):
	seqs = []
	choices = [degen_nucs[n] for n in primer]
	all_combinations_gen(0,len(primer),[],seqs,choices)
	return seqs

primers = expandprimer(primer)
print('%i primers' %len(primers))
for p in primers:
	print p

fs = Fasta.FastaSeqs()
fs.loadseqs([f])

#primers.append('A')

seqmatches = {}

for s in fs.seqs:
	seqmatches[s] = []
	for p in primers:

		prgx = re.compile(p,re.IGNORECASE)
		for match in prgx.finditer(fs.seqs[s].seq):
			seqmatches[s].append( [p,match.start(),''] )

		rvscmp = re.compile(rvs_comp_str(p),re.IGNORECASE)
		for match in rvscmp.finditer(fs.seqs[s].seq):
			seqmatches[s].append( [p,match.start(),'(rvs)'] )

for s in seqmatches:
	print(s)
	if len(seqmatches[s]) == 0:
		print('\tNO MATCHES')
		continue
	for m in seqmatches[s]:
		print('\t%s' %('\t'.join([str(i) for i in m])))
