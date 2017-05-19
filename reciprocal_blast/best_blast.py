#!/usr/bin/env python

'''
this step of reporting the best (one-directional) blast hit is brutally simple,
since blast results are already intra-query sorted and so you simply take the first one
'''

import sys

'''
Cre08.g379050.t1.1	209923	52.11	71	34	0	5	75	8	78	3e-19	90.5
Cre08.g378150.t1.1	451877	48.07	337	144	8	255	583	270	583	6e-95	  307
Cre08.g378150.t1.1	451877	45.16	93	47	2	105	197	80	168	3e-16	82.4
'''

bests = {}

for l in open(sys.argv[1]):
	f = l.strip().split()
	A = f[0]
	if not A in bests:
		print('%s\t%s' %(A,f[1]))
		bests[A] = 1
