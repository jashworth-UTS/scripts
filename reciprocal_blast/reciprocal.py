#!/usr/bin/env python
import sys

fa = sys.argv[1]
fb = sys.argv[2]

bestAs = {}
for l in open(fa):
	f = l.strip().split()
	bestAs[f[0]] = f[1]

# because reciprocal hits are required, can simply filter through 'b' now for 'a' matches
for l in open(fb):
	f = l.strip().split()
	a = f[1]
	b = f[0]
	if not a in bestAs: continue
	if bestAs[a] == b:
		print('%s\t%s' %(a,b))
