#!/usr/bin/python

import sys, string
from AshworthUtil import rvs_comp_str

if len(sys.argv) > 1: input = string.join( sys.argv[1:], '' )
else: input = sys.stdin.read()

print rvs_comp_str(input.rstrip('\n'))
