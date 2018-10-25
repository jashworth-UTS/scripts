#!/usr/bin/env python
import sys

fa = sys.argv[1]
fb = sys.argv[2]

# get top blast hits for A
bestAs = {}
for l in open(fa):
	f = l.strip().split()
	A = f[0]
	# easy because output sorted, first occurrence sid for qid (first col) is best
	if not A in bestAs:
		bestAs[A] = f[1]

# get top blast hits for B
bestBs = {}
for l in open(fb):
	f = l.strip().split()
	B = f[0]
	if not B in bestBs:
		bestBs[B] = f[1]

# output all reciprocal matches
recip = {}
for A,B in bestAs.items():
	if not B in bestBs: continue
	BA = bestBs[B]
	if not BA == A: continue
	print('%s\t%s' %(A,B))
