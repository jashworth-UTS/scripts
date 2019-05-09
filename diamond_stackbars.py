#!/usr/bin/env python3

import pandas
import matplotlib
from matplotlib import pyplot
import json
import sys

def msg(msg):
	sys.stderr.write('%s\n' %msg)

matplotlib.rcParams['hatch.linewidth'] = 0.02

# makes some nice stacked bar plots for abundance counts of annotations for diamond matches
# e.g. RNAseq reads > diamond blastx vs. uniprot_trembl > compile stats (separate script) > plot (here)

dat = json.load(open(sys.argv[1]))
stats = [
#	'os',
#	'os_gn',
#	'gn',
#	'ox',
#	'gn_acc',
	'ox_cat',
]
n_topcls = 40
normtype = 'all'
#normtype = 'incl'

group_order = []
sample_groups = {}
sg = pandas.read_csv('samples','\t')
for i,row in sg.iterrows():
	grp = row['group']
	if not grp in group_order: group_order.append(grp)
	if not grp in sample_groups: sample_groups[grp] = []
	prf = str(row['prefix'])
	sample_groups[grp].append(prf)

collapse_spc = False
if collapse_spc:
	msg('collapsing species/strains...')
	for smp in dat:
		collapsed = {}
		drop = []
		for cls,val in dat[smp]['os'].items():
			key = ' '.join(cls.strip().split()[:2])
			if not key in collapsed: collapsed[key] = 0
			collapsed[key] += val
			drop.append(cls)
		for d in drop: dat[smp]['os'].pop(d)
		for k,v in collapsed.items():
			dat[smp]['os'][k] = v

collapse_samples = True
collapse_groups = True
if collapse_groups:
	msg('collapsing groups...')
	new = {}
	for grp,prfs in sample_groups.items():
		new[grp] = {}
		for smp in dat:
			for prf in prfs:
				if smp.startswith(prf):
					for stat,kv in dat[smp].items():
						if not stat in new[grp]: new[grp][stat] = {}
						for k,v in kv.items():
							if not k in new[grp][stat]: new[grp][stat][k] = v
							else: new[grp][stat][k] += v
	dat = new
	# collapse_groups has the same end result as collapse_samples
	collapse_samples = False

if 'ox_cat' in stats:
	# combined 'NA' with 'Unclassified'
	for grp in dat:
		dat[grp]['ox_cat']['Unclassified'] += dat[grp]['ox_cat']['NA']
		dat[grp]['ox_cat'].pop('NA')

if collapse_samples:
	msg('collapsing samples...')
	new = {}
	for smp in dat:
		key = smp
		for grp,prfs in sample_groups.items():
			for prf in prfs:
				if smp.startswith(prf): key = prf
		if not key in new: new[key] = {}
		for stat,kv in dat[smp].items():
			if not stat in new[key]: new[key][stat] = {}
			for k,v in kv.items():
				if not k in new[key][stat]: new[key][stat][k] = v
				else: new[key][stat][k] += v
	dat = new

print(dat)

for stat in stats:

	msg('%s...' %stat)

	allclasses = {}
	samples = sorted(dat.keys())
	norms = {}
	maxsum = 0
	allreads = {}
	for smp in samples:
		classes = dat[smp][stat]
		n_topcls = min(n_topcls,len(classes))
		norms[smp] = {}
		valsum = 0

		if normtype == 'all':
		# normalize 'fraction total' by ALL classes
			for val in classes.values(): valsum += int(val)

		elif normtype == 'incl':
			for cls,val in sorted(classes.items(), key=lambda x: x[1], reverse=True)[:n_topcls]:
				#	normalize 'fraction total' by only included classes
				valsum += int(val)

		if valsum > maxsum: maxsum = valsum
		allreads[smp] = valsum

		for cls,val in sorted(classes.items(), key=lambda x: x[1], reverse=True)[:n_topcls]:
			norms[smp][cls] = float(val)/valsum
			if not cls in allclasses: allclasses[cls] = 0
			allclasses[cls] += int(val)

	print(allclasses)

#	plots = ['abs','rel']
#	plots = ['rel']
	plots = ['rel','vir']

	nsamples = len(dat)
	fig,axes = pyplot.subplots(len(plots),1)
	if len(plots)==1: axes = [axes]
	ax = {}
	for i in range(len(plots)): ax[plots[i]] = axes[i]
	fig.set_size_inches(8,8)
	rmarg = 0.60
	pyplot.subplots_adjust(bottom=0.15,top=0.95,hspace=0.3,left=0.15,right=rmarg)
	cmap = matplotlib.cm.get_cmap('tab10')
	ncolors = 10

	legend = []
	# legend for most abundant over all samples
	for cls,n in sorted(allclasses.items(),key=lambda x: x[1],reverse=True)[:n_topcls]:
		legend.append(cls)

	def coli_hatch(cls):
		keys = [k for k,v in sorted(allclasses.items(),key=lambda x: x[1],reverse=True)]
		#coli = float(keys.index(cls)) / len(keys)
		ki = keys.index(cls)
		coli = float(ki % ncolors)/ncolors
		hatches = ['///','\\\\\\','xxx','ooo','**','OOO','...']
		if len(allclasses) <= ncolors: hatches = ['']
		hi = round((ki - (ki % ncolors))/ncolors)
		hatch = ''
		if hi < len(hatches): hatch = hatches[hi]
		print(cls,ki,coli,hi,hatch)
		return((coli,hatch))
	
	ax[plots[0]].set_title('%i most abundant (%s)' %(n_topcls,stat))
	ax[plots[len(plots)-1]].set_xlabel('Sample')

	smporder = []
	for grp in group_order:
		for prf in sorted(sample_groups[grp]):
			if prf in dat: smporder.append(prf)
		if grp in dat: smporder.append(grp)
	xticklabs = smporder
	xticklabs_full = []
	for grp in group_order:
		for prf in sorted(sample_groups[grp]):
			if prf in dat: xticklabs_full.append('%s: %s' %(grp,prf))
		if grp in dat: xticklabs_full.append('%s' %grp)

	if 'abs' in plots:
		ax['abs'].set_ylabel('Reads')
		ax['abs'].set_ylim(0,maxsum)

		ax['abs'].set_xlim(0.5,nsamples+0.5)
		ax['abs'].set_xticks(range(1,nsamples+1))
#		ax['abs'].set_xticklabels([])
		ax['abs'].set_xticklabels(xticklabs,rotation=60,ha='right',fontsize=8)

		vals = ax['abs'].get_yticks()
		ax['abs'].set_yticklabels(['%g' %v for v in vals])

		x = 0
		for smp in smporder:
			bot = 0
			x += 1
			ax['abs'].bar(x,height=allreads[smp],bottom=bot,color='lightgray',edgecolor='black',linewidth=0.02)
			for cls,val in sorted(dat[smp][stat].items(),key=lambda x: x[1],reverse=True)[:n_topcls]:
				ch = coli_hatch(cls)
				ax['abs'].bar(x,height=val,bottom=bot,color=cmap(ch[0]),edgecolor='black',hatch=ch[1],linewidth=0.02)
				bot += val

	if 'rel' in plots:
		if not 'abs' in plots: pyplot.subplots_adjust(left=0.11)
		# norm values
		ax['rel'].set_ylabel('Proportion of Reads (within-sample)')
		ax['rel'].set_ylim(0,1)

		ax['rel'].set_xlim(0.5,nsamples+0.5)
		ax['rel'].set_xticks(range(1,nsamples+1))
		ax['rel'].set_xticklabels(xticklabs_full,rotation=60,ha='right',fontsize=8)

		thresh = 0.05
		x = 0
		for smp in smporder:
			bot = 0
			x += 1
			for cls,nrm in sorted(norms[smp].items(),key=lambda x: x[1],reverse=True)[:n_topcls]:
				ch = coli_hatch(cls)
				ax['rel'].bar(x,height=nrm,bottom=bot,color=cmap(ch[0]),edgecolor='black',hatch=ch[1],linewidth=0.02)
				bot += nrm

				# also add anything above a certain relative threshold within any sample
				if nrm>=thresh:
					if not cls in legend:
						legend.append(cls)

	if 'vir' in plots:
		# virus values
		ax['vir'].set_ylabel('Proportion of Viral Reads (within-sample)')
		cls = 'Viruses'
		ymax = 0
		for smp,cv in norms.items():
			if cv[cls]>ymax: ymax=cv[cls]
		ax['vir'].set_ylim(0,ymax)

		ax['vir'].set_xlim(0.5,nsamples+0.5)
		ax['vir'].set_xticks(range(1,nsamples+1))
		ax['vir'].set_xticklabels(xticklabs_full,rotation=60,ha='right',fontsize=8)

		x = 0
		for smp in smporder:
			bot = 0
			x += 1
			ch = coli_hatch(cls)
			ax['vir'].bar(x,height=norms[smp][cls],bottom=bot,color=cmap(ch[0]),edgecolor='black',hatch=ch[1],linewidth=0.02)

	counter = 0
	for cls in legend:
		n = allclasses[cls]
		ch = coli_hatch(cls)
		fig.text(rmarg+0.05,0.90-counter*0.02,'%s (%i)' %(cls,n),color=cmap(ch[0]),fontsize=8)
#		fig.bar(rmarg+0.02,height=0.02,bottom=0.86-counter*0.03,color=cmap(ch[0]),edgecolor='black',hatch=ch[1],linewidth=0.02)
		fig.patches.extend([
			matplotlib.patches.Rectangle((rmarg+0.01,0.895-counter*0.02),0.02,0.015,linewidth=0.02,facecolor=cmap(ch[0]),edgecolor='black',hatch=ch[1],transform=fig.transFigure,figure=fig)
		])
		counter += 1
	
	fig.savefig('diamond_stats_%s.pdf' %stat)
#	fig.savefig('diamond_stats_%s.svg' %stat)
