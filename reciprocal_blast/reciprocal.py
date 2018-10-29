#!/usr/bin/env python
import sys

fa = sys.argv[1]
fb = sys.argv[2]

taken = {}

minaln = 70
minprc = 40
maxevl = 1e-20

'''
TRINITY_DN28071_c0_g1_i1	242653	29.592	196	6.87e-22
TRINITY_DN28071_c0_g1_i1	154729	53.226	62	1.49e-18
TRINITY_DN28033_c0_g1_i1	206643	33.846	130	5.10e-20
'''

# get first(top) hits for A
bestAs = {}
for l in open(fa):
	fs = l.strip().split()
	A = fs[0]
	prc = float(fs[2])
	aln = int(fs[3])
	evl = float(fs[4])
	if prc < minprc: continue
	if aln < minaln: continue
	if evl > maxevl: continue
	# easy because blast output is sorted: first occurrence is best eval/bitscore
	if A in bestAs: continue
	bestAs[A] = fs[1]
	taken[A] = False
	taken[fs[1]] = False

# get first(top) hits for B
bestBs = {}
for l in open(fb):
	fs = l.strip().split()
	B = fs[0]
	prc = float(fs[2])
	aln = int(fs[3])
	evl = float(fs[4])
	if prc < minprc: continue
	if aln < minaln: continue
	if evl > maxevl: continue
	if B in bestBs: continue
	bestBs[B] = fs[1]
	taken[B] = False
	taken[fs[1]] = False

# output all reciprocal matches
recip = {}
for A,AB in bestAs.items():
	if taken[A]: continue
	if taken[AB]: continue
	if not AB in bestBs: continue
	if not bestBs[AB] == A: continue
	# is reciprocal best blast, mark as 'taken' and print
	taken[A] = True
	taken[AB] = True
	print('%s\t%s' %(A,AB))
