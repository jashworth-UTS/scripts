#!/usr/bin/env python

import sys


a = sys.argv[1]
b = sys.argv[2]

print('%s\t%s' %(a,b))

fa = '%s.vs.%s.best' %(a,b)
fb = '%s.vs.%s.best' %(b,a)

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
