#!/usr/bin/env python
import sys

fa = sys.argv[1]
fb = sys.argv[2]

bestAs = {}
for l in open(fa):
	f = l.strip().split()
	A = f[0]
	# avoid repeats for A
	if not A in bestAs:
		bestAs[A] = f[1]

# because reciprocal hits are required, can simply filter through 'b' now for 'a' matches
recip = {}
Btaken = []
for l in open(fb):
	f = l.strip().split()
	A = f[1]
	B = f[0]
	if B in Btaken: continue
	if not A in bestAs: continue
	if bestAs[A] == B:
		# avoid duplicate lines (such as for repeat proteins repeatedly matching to each other)
		if not A in recip:
			recip[A]=B
			Btaken.append(B)

for A,B in recip.items():	
	print('%s\t%s' %(A,B))
