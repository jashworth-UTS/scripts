#!/usr/bin/env python
import sys,re
l = filter(None, re.split("[ \t\n,]+", sys.stdin.read()))
u = []
for i in l:
	if not i in u:
		u.append(i)
print( '\n'.join([i for i in u]))
