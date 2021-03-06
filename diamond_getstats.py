#!/usr/bin/env python3

import sys
import re
import json

# OS=Pelagophyceae sp. 03-D4 OX=882321 GN=rbcL
# this regex is currently written specifically for a trembl DB
# (modify/add more regex types as needed--further code will also need modifications)
rgxs = {
	'diamond_trembl' : re.compile('tr\|([^\|]+?)\|.*?OS=(.*?) OX=(.*?) GN=([^ ]+)')
}

# only doing one regex at a time here for now
rgx = rgxs['diamond_trembl']

fls = sys.argv[1:]

max_incl = 1

stats = {}
for fl in fls:
	stats[fl] = {
		'os':{},
#		'ox':{},
		'gn':{},
		'gn_acc':{},
	}
	counts = {}
	for l in open(fl):
		qid = l.split()[0]

		# limit included tabulation to first (top) N matches
		if not qid in counts: counts[qid] = 0
		counts[qid] += 1
		if counts[qid] > max_incl: continue

		m = rgx.search(l)
		if m == None: continue
		acc,os,ox,gn = m.groups()
		gn_acc = '%s_%s' %(gn,acc)
		if not os in stats[fl]['os']: stats[fl]['os'][os] = 0
#		if not ox in stats[fl]['ox']: stats[fl]['ox'][ox] = 0
		if not gn in stats[fl]['gn']: stats[fl]['gn'][gn] = 0
		if not gn_acc in stats[fl]['gn_acc']: stats[fl]['gn_acc'][gn_acc] = 0
		stats[fl]['os'][os] += 1
#		stats[fl]['ox'][ox] += 1
		stats[fl]['gn'][gn] += 1
		stats[fl]['gn_acc'][gn_acc] += 1

json.dump(stats,open('diamond_stats.json','w'))
