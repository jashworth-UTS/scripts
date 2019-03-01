#!/usr/bin/env python3

# plots some pre-parsed PBSPro qstat -f stats from a custom logger
# for visualizing cpu and mem usage over time

import sys
import re
import matplotlib
from matplotlib import pyplot

headerfile = sys.argv[1]
statsfile = sys.argv[2]

delim = '\t'

fields = open(headerfile).read().strip().split(delim)

stats = {}
mintime = None
for l in open(statsfile):
	ll = l.strip().split(delim)
	jobid = ll[fields.index('JobId')]
	time = ll[fields.index('time')]
	if mintime == None or time < mintime: mintime = time
#	time = (ll[fields.index('time')], ll[fields.index('resources_used.walltime')])
#	if mintime == None or time[0] < mintime: mintime = time[0]
	if not jobid in stats: stats[jobid] = {}
	if not time in stats[jobid]: stats[jobid][time] = {}
	for i in range(len(fields)):
		fld = fields[i]
		val = ll[i]
		if 'kb' in val: val = float(re.sub('kb','',val))*1e3
		elif 'mb' in val: val = float(re.sub('mb','',val))*1e6
		elif 'gb' in val: val = float(re.sub('gb','',val))*1e9
		stats[jobid][time][fld] = val

plotfields = [
#	'Resource_List.mem',
#	'Resource_List.ncpus',
	'resources_used.cpupercent',
#	'resources_used.cput',
	'resources_used.mem',
#	'resources_used.ncpus',
#	'resources_used.vmem',
]

cmap = matplotlib.cm.get_cmap('tab10')
ncolors = 10

fig,axes = pyplot.subplots(len(plotfields),1)
#fig.set_size_inches(6,6)
rmarg = 0.7
pyplot.subplots_adjust(hspace=0.3,bottom=0.1,top=0.95,right=rmarg)

jobids = sorted(stats)
def coli(jobid):
	return(float(jobids.index(jobid) % ncolors)/ncolors)

# currently system time is in seconds and this converts to minutes
def tfmt(t):
	return(float(int(t)-int(mintime))/60)

for i in range(len(plotfields)):
	fld = plotfields[i]
	maxy = 0
	for jobid,job in sorted(stats.items(), key=lambda x: x[0]):
		x = []
		y = []
		for t in sorted(stats[jobid], key=lambda x: int(x[0])):
			x.append(tfmt(t))
			y.append(float(stats[jobid][t][fld]))
		axes[i].plot(x,y,label=jobid,color=cmap(coli(jobid)))
		maxy = max(maxy,max(y))

		# add a dashed line for the requested/allocated memory for this job
		if fld in ['resources_used.mem','resources_used.vmem']:
			x = []
			y = []
			for t in sorted(stats[jobid], key=lambda x: int(x[0])):
				x.append(tfmt(t))
				y.append(float(stats[jobid][t]['Resource_List.mem']))
			axes[i].plot(x,y,label=jobid,color=cmap(coli(jobid)),dashes=(2,2))
			maxy = max(maxy,max(y))

		# add a dashed line for the requested/allocated ncpus for this job
		if fld in ['resources_used.ncpus']:
			x = []
			y = []
			for t in sorted(stats[jobid], key=lambda x: int(x[0])):
				x.append(tfmt(t))
				y.append(float(stats[jobid][t]['Resource_List.ncpus']))
			axes[i].plot(x,y,label=jobid,color=cmap(coli(jobid)),dashes=(2,2))
			maxy = max(maxy,max(y))

		# add a dashed line for the max 'cpupercent' available for the requested/allocated ncpus
		if fld in ['resources_used.cpupercent']:
			x = []
			y = []
			for t in sorted(stats[jobid], key=lambda x: int(x[0])):
				x.append(tfmt(t))
				y.append(float(stats[jobid][t]['Resource_List.ncpus'])*100)
			axes[i].plot(x,y,label=jobid,color=cmap(coli(jobid)),dashes=(2,2))
			maxy = max(maxy,max(y))

		axes[i].set_ylim(0,maxy*1.1)
		vals = axes[i].get_yticks()
		axes[i].set_yticklabels(['%g' %v for v in vals])

	axes[i].set_ylabel(fld)
	axes[i].set_xlabel('time (min)')

for i in range(len(jobids)):
	fig.text(rmarg+0.05,0.9-i*0.05, jobids[i], color=cmap(coli(jobids[i])))
	
fig.savefig('stats.pdf')
