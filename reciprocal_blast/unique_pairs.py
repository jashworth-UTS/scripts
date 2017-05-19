#!/usr/bin/env python
import sys,re
fs = sys.argv[1:]
for i in range(len(fs)):
	ni = fs[i]
	for j in range(i+1,len(fs)):
		nj = fs[j]
		ni = re.sub('.aa.ids.fa','',ni)
		nj = re.sub('.aa.ids.fa','',nj)
		print('%s\t%s' %(ni,nj))
