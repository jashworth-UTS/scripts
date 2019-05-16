#!/usr/bin/env python3

import sys
import re
from Fasta import FastaSeqs
from AshworthUtil import rvs_comp_str

# braker didn't output spliced CDS sequences, so parse the gtf file and create them

def msg(msg): sys.stderr.write('%s\n' %msg)

class Gene(dict):
	def __init__(self):
		dict.__init__(self)

parts = [
	'CDS',
]

re_parent = re.compile('transcript_id ([^;]+)')

genes = {}

gtf = sys.argv[1]
for l in open(gtf):
	ll = l.strip().split('\t')
	type = ll[2]
	if type in parts:
		seq = ll[0]
		start = int(ll[3])
		end = int(ll[4])
		strand = ll[6]
		parent = re_parent.search(ll[8]).groups()[0]
		if not parent in genes:
			genes[parent] = Gene()
			genes[parent]['seq'] = seq
			genes[parent]['cds'] = []
		genes[parent]['cds'].append( (start,end,strand) )

seqfile = sys.argv[2]
src = FastaSeqs()
src.loadseqs([seqfile])

outp = []

for id,gene in genes.items():
	seq = gene['seq']
	if not seq in src.seqs:
		msg('%s not found in source!' %seq)
		sys.exit()
	if 'cds' in gene:
		fullcds = []
		rvs = False
		for start,end,strand in sorted(gene['cds'],key=lambda x: x[0]):
			ss = src.seqs[seq].seq[(start-1):end]
			if strand == '-':
				ss = rvs_comp_str(ss)
				rvs = True
			fullcds.append(ss)
		if rvs: fullcds.reverse()
		outp.append('>%s_%s\n%s' %(id,'cds',''.join(fullcds)))

print('\n'.join(outp))
