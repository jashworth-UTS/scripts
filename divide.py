#!/usr/bin/env python

import sys
import Fasta

fs = sys.argv[1:]

npf = 1000

for f in fs:
	ss = Fasta.FastaSeqs()
	ss.loadseqs([f])
	ns = len(ss.seqs)
	for i in range(0,ns,npf):
		of = open('%s.%i.fa' %(f,i),'w')
		end = min(i+npf,ns)
		keys = ss.order[i:end]
		for key in keys:
			of.write('%s\n' %ss.seqs[key].fasta())
		of.close()
