# J Ashworth 2020 April

# Plot comparative distributions of expression data for selected transcripts

import numpy
import sys
import random

import pandas as pd
#tb = pd.read_table("Phatr.data.test", header=0, index_col=0, sep=" ")
#tb = pd.read_table("Phatr.data", header=0, index_col=0, sep=" ")
tb = pd.read_table("Phatr.data.RNAseq", header=0, index_col=0, sep=" ")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
fig,ax = plt.subplots()

txs = [
	"Phatr3_J54958",
	"Phatr3_J45582",
	"Phatr3_J14260",
	"Phatr3_J16798",
	"Phatr3_J35766",
	"Phatr3_J49202",
	"Phatr3_J18049",
]

#overall median (all transcripts)
mdall = numpy.nanmedian(tb.values)
ax.add_line(Line2D([-1,len(txs)],[mdall,mdall],linestyle='--',c='black'))

ax.set_xlim(-0.5,len(txs)-0.5)
ax.set_xticks(list(range(len(txs))))
ax.set_ylabel("Normalized Expression (TPM)")
spread = 0.2

# get stats first, for sorted output
ranks = []
for tx in txs:
	ranks.append( (tx, numpy.nanmedian(tb.loc[tx]), numpy.nanstd(tb.loc[tx])) )

xticklabs = []
i = 0
for tx,md,sd in sorted(ranks, key=lambda r: r[2], reverse=True):
	xticklabs.append(tx)
	for j in range(len(tb.columns)):
		ax.scatter(i+random.uniform(-1*spread,spread),tb.loc[tx][j],c="black",s=5,alpha=0.3,edgecolor='none')
	sys.stderr.write('%s %0.2f\n' %(tx,md))
	ax.add_line(Line2D([i-spread-0.1,i+spread+0.1],[md,md],c="black",linewidth=2))
	i += 1

ax.set_xticklabels(xticklabs,rotation=45,ha="right")

plt.subplots_adjust(0.1,0.2,0.95,0.95)

fig.savefig('fig.pdf')
#fig.savefig('fig.png')
plt.close()
