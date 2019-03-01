#!/usr/bin/env python3

import matplotlib
from matplotlib import pyplot
import json
import sys

# makes some nice stacked bar plots for abundance counts of annotations for diamond matches
# e.g. RNAseq reads > diamond blastx vs. uniprot_trembl > compile stats (separate script) > plot (here)

dat = json.load(open(sys.argv[1]))
fields = [
	'os',
	'gn',
	'gn_acc',
]
n_topcls = 30

for field in fields:

	allclasses = {}
	grps = sorted(dat.keys())
	norms = {}
	maxsum = 0
	for grp in grps:
		classes = dat[grp][field]
		norms[grp] = {}
		valsum = 0
		for val in classes.values(): valsum += val
		if valsum > maxsum: maxsum = valsum
		for cls,val in sorted(classes.items(), key=lambda x: x[1], reverse=True)[:n_topcls]:
			norms[grp][cls] = float(val)/valsum
			if not cls in allclasses: allclasses[cls] = 0
			allclasses[cls] += val

	ngrps = len(dat)
	fig,axes = pyplot.subplots(2,1)
	axes[0].set_title('%i most abundant (%s)' %(n_topcls,field))
	rmarg = 0.5
	pyplot.subplots_adjust(bottom=0.3,left=0.15,right=rmarg)
	cmap = matplotlib.cm.get_cmap('tab10')

	# abs values
	axes[0].set_ylabel('Reads')
	axes[0].set_ylim(0,maxsum)
	axes[0].set_xlim(0.5,ngrps+0.5)
	axes[0].set_xticks(range(1,ngrps+1))
	axes[0].set_xticklabels([])

	def coli(cls):
		keys = [k for k,v in sorted(allclasses.items(),key=lambda x: x[1])]
		#coli = float(keys.index(cls)) / len(keys)
		coli = float(keys.index(cls) % 10)/10
		print(cls,coli)
		return(coli)

	x = 0
	for grp in grps:
		bot = 0
		x += 1
		for cls,val in sorted(dat[grp][field].items(),key=lambda x: x[1],reverse=True)[:n_topcls]:
			axes[0].bar(x,height=val,bottom=bot,color=cmap(coli(cls)))
			bot += val

	# norm values
	axes[1].set_ylabel('Fraction of Reads')
	axes[1].set_ylim(0,1)
	axes[1].set_xlim(0.5,ngrps+0.5)
	axes[1].set_xticks(range(1,ngrps+1))
	xticklabs = [g[:10]+'..'+g[-3:] for g in grps]
	axes[1].set_xticklabels(xticklabs,rotation=60,ha='right')

	x = 0
	for grp in grps:
		bot = 0
		x += 1
		for cls,nrm in sorted(norms[grp].items(),key=lambda x: x[1],reverse=True)[:n_topcls]:
			axes[1].bar(x,height=nrm,bottom=bot,color=cmap(coli(cls)))
			bot += nrm

	counter = 0
	for cls,n in sorted(allclasses.items(),key=lambda x: x[1],reverse=True)[:n_topcls]:
		fig.text(rmarg+0.05,0.9-counter*0.03,'%s (%i)' %(cls,n),color=cmap(coli(cls)))
		counter += 1

	fig.savefig('diamond_stats_%s.pdf' %field)
