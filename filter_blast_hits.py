#!/usr/bin/env python

import sys
import pandas as pd

blastfile = sys.argv[1]
tb = pd.read_csv(blastfile,sep='\t',comment='#',index_col=False,header=None)
tb = tb[ (tb[2] > 35) & (tb[3] > 50) ]
ids = tb[1]
print('\n'.join(ids))
