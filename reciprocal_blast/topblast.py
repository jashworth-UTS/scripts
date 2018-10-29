#!/usr/bin/env python

import sys

'''
TRINITY_DN28071_c0_g1_i1	38647	33.019	212	6.57e-29
TRINITY_DN28012_c0_g1_i1	9790	99.438	178	2.54e-123
TRINITY_DN28055_c0_g1_i1	5273	100.000	69	6.98e-42
TRINITY_DN28055_c0_g2_i1	5273	100.000	83	1.40e-51
'''

qs = {}

for l in open(sys.argv[1]):
	fs = l.strip().split()
	qid = fs[0]
	sid = fs[1]
	prc = float(fs[2])
	aln = int(fs[3])
	evl = float(fs[4])
	if qid in qs: continue # only take first hit encountered since blast results are sorted
	qs[qid] = (sid, prc, aln, evl)

sys.stderr.write('%i top hits\n' %len(qs))

for q,s in qs.items():
	print('\t'.join([s[0],str(int(s[1])),str(s[2]),'%0.2g' %s[3]]))
