#!/usr/bin/env python3

import matplotlib
from matplotlib import pyplot
import json
import sys

matplotlib.rcParams['hatch.linewidth'] = 0.02

# makes some nice stacked bar plots for abundance counts of annotations for diamond matches
# e.g. RNAseq reads > diamond blastx vs. uniprot_trembl > compile stats (separate script) > plot (here)

dat = json.load(open(sys.argv[1]))
fields = [
	'os',
	'gn',
	'ox',
#	'gn_acc',
]
n_topcls = 45
normtype = 'all'
#normtype = 'incl'

collapse_spc = True

if collapse_spc:
	for grp in dat:
		collapsed = {}
		drop = []
		for cls,val in dat[grp]['os'].items():
			key = ' '.join(cls.strip().split()[:2])
			if not key in collapsed: collapsed[key] = 0
			collapsed[key] += val
			drop.append(cls)
		for d in drop: dat[grp]['os'].pop(d)
		for k,v in collapsed.items():
			dat[grp]['os'][k] = v

for field in fields:

	allclasses = {}
	grps = sorted(dat.keys())
	norms = {}
	maxsum = 0
	for grp in grps:
		classes = dat[grp][field]
		norms[grp] = {}
		valsum = 0

		if normtype == 'all':
		# normalize 'fraction total' by ALL classes
			for val in classes.values(): valsum += int(val)
		elif normtype == 'incl':
			for cls,val in sorted(classes.items(), key=lambda x: x[1], reverse=True)[:n_topcls]:
				#	normalize 'fraction total' by only included classes
				valsum += int(val)
		if valsum > maxsum: maxsum = valsum

		for cls,val in sorted(classes.items(), key=lambda x: x[1], reverse=True)[:n_topcls]:
			norms[grp][cls] = float(val)/valsum
			if not cls in allclasses: allclasses[cls] = 0
			allclasses[cls] += int(val)

	ngrps = len(dat)
	fig,axes = pyplot.subplots(2,1)
	fig.set_size_inches(8,8)
	axes[0].set_title('%i most abundant (%s)' %(n_topcls,field))
	rmarg = 0.65
	pyplot.subplots_adjust(bottom=0.2,left=0.15,right=rmarg)
	cmap = matplotlib.cm.get_cmap('tab10')
	ncolors = 10

	# abs values
	axes[0].set_ylabel('Reads')
	axes[0].set_ylim(0,maxsum)
	axes[0].set_xlim(0.5,ngrps+0.5)
	axes[0].set_xticks(range(1,ngrps+1))
	axes[0].set_xticklabels([])

	def coli_hatch(cls):
		keys = [k for k,v in sorted(allclasses.items(),key=lambda x: x[1],reverse=True)]
		#coli = float(keys.index(cls)) / len(keys)
		ki = keys.index(cls)
		coli = float(ki % ncolors)/ncolors
		hatches = ['///','\\\\\\','xxx','ooo','***','OOO','...']
		hi = round((ki - (ki % ncolors))/ncolors)
		hatch = ''
		if hi < len(hatches): hatch = hatches[hi]
		print(cls,ki,coli,hi,hatch)
		return((coli,hatch))
	
	x = 0
	for grp in grps:
		bot = 0
		x += 1
		for cls,val in sorted(dat[grp][field].items(),key=lambda x: x[1],reverse=True)[:n_topcls]:
			ch = coli_hatch(cls)
			axes[0].bar(x,height=val,bottom=bot,color=cmap(ch[0]),edgecolor='black',hatch=ch[1],linewidth=0.02)
			bot += val

	# norm values
	axes[1].set_ylabel('Fraction of Reads')
	axes[1].set_ylim(0,1)
	axes[1].set_xlim(0.5,ngrps+0.5)
	axes[1].set_xticks(range(1,ngrps+1))
	xticklabs = [g[:15] for g in grps]
#	xticklabs = [g[:10]+'..'+g[-3:] for g in grps]
	axes[1].set_xticklabels(xticklabs,rotation=60,ha='right',fontsize=8)

	x = 0
	for grp in grps:
		bot = 0
		x += 1
		for cls,nrm in sorted(norms[grp].items(),key=lambda x: x[1],reverse=True)[:n_topcls]:
			ch = coli_hatch(cls)
			axes[1].bar(x,height=nrm,bottom=bot,color=cmap(ch[0]),edgecolor='black',hatch=ch[1],linewidth=0.02)
			bot += nrm

	counter = 0
	for cls,n in sorted(allclasses.items(),key=lambda x: x[1],reverse=True)[:n_topcls]:
		ch = coli_hatch(cls)
		fig.text(rmarg+0.05,0.9-counter*0.02,'%s (%i)' %(cls,n),color=cmap(ch[0]),fontsize=8)
#		fig.bar(rmarg+0.02,height=0.02,bottom=0.86-counter*0.03,color=cmap(ch[0]),edgecolor='black',hatch=ch[1],linewidth=0.02)
		fig.patches.extend([
			matplotlib.patches.Rectangle((rmarg+0.01,0.895-counter*0.02),0.02,0.015,linewidth=0.02,facecolor=cmap(ch[0]),edgecolor='black',hatch=ch[1],transform=fig.transFigure,figure=fig)
		])
		counter += 1

	fig.savefig('diamond_stats_%s.pdf' %field)
