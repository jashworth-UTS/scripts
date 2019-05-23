#!/usr/bin/env python3

'''
Python utility script for reading and auto-processing .fcs flow cytometry data files.
See required dependencies (library imports) below: install them via 'sudo pip3 install [x]' under linux or OS X
# (not tested on Windows)
J Ashworth May 2019
Work in progress: fork/expand/improve as appropriate
'''

import FlowCytometryTools
from FlowCytometryTools import FCMeasurement
import matplotlib
from matplotlib import pyplot as plt
cmap = matplotlib.cm.get_cmap('tab10')
ncols = 10
def cmap_col(i): return(cmap(float(i % ncols)/ncols))
from sys import stderr,argv
import scipy.stats
import numpy
from sklearn.cluster import *
import re

def msg(msg): stderr.write('%s\n' %msg)

class Cyto:
	def __init__(self,fcs):
		self.prf = fcs
		msg('loading %s' %fcs)
		self.dat = FCMeasurement(ID=fcs,datafile=fcs)
		self.clusters = {}
		self.stats = {}

	def cluster(self,clusterer,chdefs=[]):
		chs = []
		for ch,scale in chdefs: chs.append(ch)
		cldat = self.dat.data[chs]
		for ch,scale in chdefs: 
			if scale == 'log': 
				# clusterers can choke on 'nan' values
				cldat[ch] = [numpy.log10(d) if d>0 else 0 for d in cldat[ch]]
		cls = clusterer.fit(cldat)

		# re-rank the clusters according to median values in first clustering dimension
		clis = set(cls.labels_)
		msg('clustering produced %i clusters:\n%s' %(len(clis),','.join([str(i) for i in clis])))
		rank = []
		for cli in sorted(clis):
#			if cli < 0: continue
			ind = [i for i in map(lambda x: x == cli, cls.labels_)]
			sortval = numpy.median(cldat[chs[0]][ind])
			rank.append( (cli,sortval,ind) )
		rank.sort(key=lambda x: x[1])
		
		for clrank in range(len(rank)):
			cli,sortval,ind = rank[clrank]
			msg('cluster %i (re-ranking to %i): size %i sortval %g' %(cli,clrank,ind.count(True),sortval))
			self.clusters['cl_%i' %clrank] = ind

	def fit_FSC_filter(self,frc=0.3,nse=12,hch='FSC-H',ach='FSC-A'):
		# problem: clumps and detritus register anomalously high FSC A/H ratios
		# sometimes, this is present mostly in high H events (but not always)
		# one approach is to try a linear fit to the lower intensity proportion of the population
		msg('fitting FSC H/A...')

		# solution 1: fit a line to the lower FSC H(/A) data; exclude any high FSC H outside of this linear relationship
		hmin = self.dat.data[hch].min()
		hmax = self.dat.data[hch].max()
		hcut = hmin + frc*(hmax-hmin)

#		# solution 1(b): as a proportion of events
#		h_ord = sorted(self.dat.data[hch])
#		h_cut = h_ord[round(float(len(h_ord))*0.75)]

		msg('h max %g h cut %g' %(hmax, hcut))
		h_fit = self.dat.data[hch][ self.dat.data[hch] <= hcut ]
		a_fit = self.dat.data[ach][ self.dat.data[hch] <= hcut ]
		lf = scipy.stats.linregress(a_fit,h_fit)

		fg,ax = plt.subplots()
		plt.subplots_adjust(left=0.2)
		ax.set_xscale('symlog')
		ax.set_yscale('symlog')
		ax.scatter(self.dat.data[ach],self.dat.data[hch],s=1,c='black',alpha=0.4)
		xv = numpy.linspace(self.dat.data[ach].min(),self.dat.data[ach].max(),50)
		ax.plot(xv,[hcut for x in xv],color='black')
		ax.set_xlabel(ach)
		ax.set_ylabel(hch)
		yv = [lf.slope*x +lf.intercept for x in xv]
		ax.plot(xv,yv,color='black',dashes=[2,2])
		se = lf.stderr
		ax.plot(xv,[(lf.slope-nse*se)*x +lf.intercept for x in xv],color='red')
#		ax.plot(xv,[(lf.slope+nse*se)*x +lf.intercept for x in xv],color='red')
		ax.set_ylim(hmin,hmax)

		fg.savefig('%s.FSCfit.png' %(self.prf))

		def se_filter(row):
			# is the height smaller than expected for this area (clumping/aggregates/detritus?)
			if row[hch] < row[ach] * (lf.slope-se*nse) + lf.intercept: return(True)
			# is the height higher than expected ('sharper'?)
#			if row[hch] > row[ach] * (lf.slope+se*nse) + lf.intercept: return(True)
			return(False)
		outliers = self.dat.data.apply(se_filter, axis=1)
		self.clusters['FSC_lf_outliers'] = outliers
		msg('%i events flagged as outliers by FSC H/A filter' %list(outliers).count(True))

	def plothist(self,ch,xrng,logdata=True,plot=True):
		if plot:
			fg,ax = plt.subplots()
			plt.subplots_adjust(left=0.2)

		clcnt = 0
		maxn = 0
		nall = self.dat.data[ch].shape[0]
		for cl in sorted(self.clusters):
			dd = self.dat.data[ch][self.clusters[cl]]
			cl_col=cmap_col(clcnt)
			clcnt += 1
			rawmed = numpy.median(dd)
			statkey = 'median_%s_%s' %(ch,cl)
			self.stats[statkey] = rawmed

			if plot:
				med = rawmed
				if logdata:
					dd = [numpy.log10(d) for d in dd if d>0]
					med = numpy.log10(med)
				nn,bins,patches = ax.hist(dd,bins=100,color=cl_col,alpha=0.3)
				n = numpy.sum(nn)
				frac = float(n)/nall
				nnmax = numpy.max(nn)
				maxn = numpy.max((maxn,nnmax))

				ax.plot([med,med],[0,1e4],color=cl_col,label='%0.2g (%0.2g)' %(rawmed,frac))

		if plot:
			ax.set_ylim(0,maxn)
			ax.legend()
			ax.set_xlabel(ch)
			ax.set_ylabel('counts')
			if not xrng == []:
				ax.set_xlim(numpy.log10(xrng[0]),numpy.log10(xrng[1]))
			fg.savefig('%s.%s.hist.png' %(self.prf,ch))

	def plot2ch(self,chs,xrng,yrng):
		msg('plotting %s' %','.join(chs))
		fg,ax = plt.subplots()
		plt.subplots_adjust(left=0.2)

		if len(self.clusters) == 0:
			ax.scatter(self.dat.data[chs[0]],self.dat.data[chs[1]],s=1,c='black',alpha=0.3)

		cls = self.clusters.keys()
#		cls = [
#			('FSC_lf_outliers','red'),
#		]
			
		clcnt = 0
		for cl in sorted(cls):
			if not cl in self.clusters: continue
			clx = self.dat.data[chs[0]][self.clusters[cl]]
			cly = self.dat.data[chs[1]][self.clusters[cl]]
			ax.scatter(clx,cly,s=1,c=cmap_col(clcnt),alpha=0.3)
			clcnt += 1

		exclude_for_ranges = 'FSC_lf_outliers'
		if exclude_for_ranges in self.clusters:
			goodx = self.dat.data[chs[0]][~self.clusters[exclude_for_ranges]]
			goody = self.dat.data[chs[1]][~self.clusters[exclude_for_ranges]]
			ax.set_xlim(goodx.min(), goodx.max())
			ax.set_ylim(goody.min(), goody.max())

		if not xrng == []: ax.set_xlim(xrng[0],xrng[1])
		if not yrng == []: ax.set_ylim(yrng[0],yrng[1])

		ax.set_xlabel(chs[0])
		ax.set_ylabel(chs[1])

#		linaxes = ['FSC','SSC','FITC','phyll']
#		linaxes = ['FSC','SSC']
		linaxes = []
		scaling = 'symlog'
		for rgx in linaxes:
			if rgx in chs[0]:
				scaling = 'linear'
				break
		ax.set_xscale(scaling)

		scaling = 'symlog'
		for rgx in linaxes:
			if rgx in chs[1]:
				scaling = 'linear'
				break
		ax.set_yscale(scaling)

		fg.savefig('%s.%s.%s.png' %(self.prf,chs[0],chs[1]))

	def allch(self,suff=''):
		msg('plotting all (histograms)')
		chs = self.dat.data.columns
		fg,axs = plt.subplots(nrows=len(chs))
		plt.subplots_adjust(left=0.2)
		for i in range(len(chs)):
			axs[i].hist(self.dat.data[chs[i]],bins=100,color='black')
			axs[i].set_ylabel(chs[i],rotation=0,ha='right')
			axs[i].set_yticklabels([])
			axs[i].set_xticklabels([])
		fg.savefig('%s.all%s.png' %(self.prf,suff))
	
	def chs(self): return(self.dat.data.columns)

	## code above is importable in ipython, jupyter, etc ##

	## code below runs if this script is called on command line e.g.:
	## `./Cyto.py *fcs`
if __name__ == "__main__":

	plots2d = [
#		(['FSC-A','SSC-A'],[0,1e7],[0,1.5e6]),
#		(['SSC-A','SSC-H'],[0,1.5e6],[0,7e5]),
#		(['FSC-A','FSC-H'],[0,1e7],[]),
		(['FSC-A','SSC-A'],[],[]),
#		(['SSC-A','SSC-H'],[],[]),
		(['FSC-A','FSC-H'],[],[]),
#		(['FSC-H','SSC-H'],[],[]),
		(['mVenus FITC-A','Chorophyll PC5.5-A'],[1e2,1e6],[1e2,1e6]),
#		(['mVenus FITC-H','Chorophyll PC5.5-H'],[1e2,1e6],[1e2,1e6]),
#		(['mVenus FITC-A','Chorophyll PC5.5-A'],[],[]),
#		(['mVenus FITC-H','Chorophyll PC5.5-H'],[],[]),
#		(['mVenus FITC-A','PE-A'],[1e2,1e7],[]),
	]

	hists = [
		('mVenus FITC-A',[1e2,1e6]),
		('Chorophyll PC5.5-A',[1e2,1e6]),
#		('mVenus FITC-A',[]),
#		('Chorophyll PC5.5-A',[]),
	]

	fs = argv[1:]
	cytos = {}
	min_events = 100
	for f in fs:
		c = Cyto(f)
		print(c.chs())

		if c.dat.data.shape[0] < min_events: continue

#		# plot 1-D histograms of all channels/data
#		c.allch()

#		# apply a linear fit filter for anomalously high FSC-A/FSC-H events
#		c.fit_FSC_filter()

		cluster = True
		if cluster:
			# use a clustering class from sklearn.cluster (e.g. KMeans) to find clusters in [custom] multi-dimensional space
			# under the assumption that one or more contains distinct, anomalous or undesireable events
			cl_chs = [
				('mVenus FITC-A','log'),
				('mVenus FITC-H','log'),
				('Chorophyll PC5.5-A','log'),
				('Chorophyll PC5.5-H','log'),
				('SSC-A','log'),
				('SSC-H','log'),
				('FSC-A','log'),
#				('FSC-H','log'),
				('Violet SSC-A','log'),
				('PE-A','log'),
#				('Time','linear'),
			]
			nclust = 5
			clusterer = KMeans(n_clusters=nclust) # simple but doesn't really find "patterns" as much as "regions"
#			clusterer = DBSCAN(eps=0.5,min_samples=100) # touchy
#			clusterer = AffinityPropagation() # too slow
#			clusterer = AgglomerativeClustering(n_clusters=nclust) # slow-ish
#			clusterer = Birch(n_clusters=nclust)
#			clusterer = MiniBatchKMeans(n_clusters=nclust)
#			clusterer = MeanShift() # too slow
#			clusterer = SpectralClustering(n_clusters=nclust) # too slow

			c.cluster(clusterer,cl_chs)
		for chs,xrng,yrng in plots2d: c.plot2ch(chs,xrng,yrng)
		for ch,xrng in hists: c.plothist(ch,xrng,plot=True)
		cytos[f] = c

	stats = []
	statkeys = []
	for c in cytos.values(): statkeys += c.stats.keys()
	statkeys = sorted(set(statkeys))
	header = ['name'] + statkeys
	print('\t'.join(header))
	res = []
	for f,c in cytos.items():
		row = [f]
		for stat in statkeys:
			val = float('nan')
			if stat in c.stats: val = c.stats[stat]
			row.append(val)
		res.append(row)
		print('\t'.join([row[0]] + ['%g' %v for v in row[1:]]))

	statplots2d = [
#		('median_mVenus FITC-A','median_Chorophyll PC5.5-A',[1e2,1e6],[1e2,1e6]),
		('median_mVenus FITC-A','median_Chorophyll PC5.5-A',[],[]),
	]
	
	namecols = {
#		'A' : 'red',
#		'B' : 'red',
#		'C' : 'red',
#		'D' : 'red',
#		'E' : 'blue',
#		'F' : 'blue',
#		'G' : 'blue',
#		'H' : 'blue',
#		'A' : 'red',
#		'B' : 'green',
#		'C' : 'blue',
#		'D' : 'gray',
#		'E' : 'purple',
#		'F' : 'orange',
#		'G' : 'cyan',
#		'H' : 'magenta',
		'A[123456]' : 'red',
		'B[123456]' : 'blue',
		'C[123456]' : 'green',
		'D[123456]' : 'orange',
	}

	for xch,ych,xrng,yrng in statplots2d:
		lab = []
		cols = []
		xv = []
		yv = []
		xi = header.index(xch)
		yi = header.index(ych)
		for row in res:
			lb = re.sub('\.fcs.*','',row[0])
			lb = lb.split('-')[-1]
			col = 'black'
			for rgx in namecols:
				if re.search(rgx,lb): col = namecols[rgx]
			cols.append(col)
			lab.append(lb)
			xv.append(row[xi])
			yv.append(row[yi])
		print(xv)
		print(yv)
		print(lab)
		fg,ax = plt.subplots()
		plt.subplots_adjust(left=0.3)
		ax.scatter(xv,yv,color=cols)
		for i in range(len(xv)):
			ax.text(xv[i],yv[i],lab[i],ha='right',color=cols[i])
		ax.set_xlabel(xch)
		ax.set_ylabel(ych)
		ax.set_xscale('log')
		ax.set_yscale('log')
		if not xrng == []: ax.set_xlim(xrng[0],xrng[1])
		if not yrng == []: ax.set_ylim(yrng[0],yrng[1])

		fg.savefig('statplot.png')
