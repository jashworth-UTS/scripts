#!/usr/bin/env python3

# J Ashworth 2019

# queries and parses PBSPro qstat -f info for simplied and/or flat job stats, helpful resource usage plots

import time
import json
import subprocess
import re

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot

# inherits from dict(), indexed by time Job[t] = {[stats at time t]}
class Job(dict):
	def __init__(self,jobid):
		dict.__init__(self)
		self.jobid = jobid

	def flatstr(self,fields,delim='\t'):
		l = []
		for t,stats in self.items():
			ll = [self.jobid,str(t)]
			for k in fields:
				v = stats[k]
				if isinstance(v,float): v = '%g' %v
				if isinstance(v,int): v = '%i' %v
				ll.append(v)
			l.append(delim.join(ll))
		return('\n'.join(l))

# runs indefinitely [e.g. in background], logging qstat statistics for jobs to file(s)
class JobLogger:

	def __init__(self,kv_delim=' = '):
		self.kv_delim = kv_delim
		self.jobs = {}

		self.fields = [
		#	'Job_Name',
		#	'Job_Owner',
			'resources_used.cpupercent',
			'resources_used.cput',
			'resources_used.mem',
			'resources_used.ncpus',
			'resources_used.vmem',
			'resources_used.walltime',
		#	'job_state',
		#	'queue',
		#	'server',
		#	'Checkpoint',
		#	'ctime',
		#	'Error_Path',
		#	'exec_host',
		#	'exec_vnode',
		#	'Hold_Types',
		#	'Join_Path',
		#	'Keep_Files',
		#	'Mail_Points',
		#	'Mail_Users',
		#	'mtime',
		#	'Output_Path',
		#	'Priority',
		#	'qtime',
		#	'Rerunable',
			'Resource_List.mem',
			'Resource_List.ncpus',
		#	'Resource_List.nodect',
		#	'Resource_List.place',
		#	'Resource_List.select',
		#	'Resource_List.walltime',
		#	'stime',
		#	'comment',
		]

	def log(self,dt,maxtime=1e6):
		self.jobs = {}
		self.tstart = round(time.time())
		logtime = 0
		while logtime < maxtime:
			logtime = time.time() - self.tstart
			self.getstats()
			self.writeflat()
			self.plot()
			time.sleep(dt)

		errf = open('jobstats.err','a')
		errf.write('max time %i elasped for logging qstat\n' %maxtime)
		errf.close()

	def getstats(self):
		qstat = subprocess.Popen(['qstat -f'],shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		ret = qstat.wait()
		stdout,stderr = qstat.communicate()
#		print(stdout,stderr) # debug
		if ret != 0:
			errf = open('jobstats.err','a')
			errf.write('%s\n' %stderr)
			errf.close()
			return()

		jobid = None
		t = round(time.time())
		for l in stdout.decode('utf-8').split('\n'):
			if l.startswith('Job Id'):
				jobid = l.strip().split()[2]
				if not jobid in self.jobs: self.jobs[jobid] = Job(jobid)
				if not t in self.jobs[jobid]: self.jobs[jobid][t] = {}
			elif self.kv_delim in l:
				kv = l.strip().split(self.kv_delim)
				if kv[0] in self.fields:
					val = kv[1]
					if 'kb' in val: val = float(re.sub('kb','',val))*1e3
					elif 'mb' in val: val = float(re.sub('mb','',val))*1e6
					elif 'gb' in val: val = float(re.sub('gb','',val))*1e9
					self.jobs[jobid][t][kv[0]] = val

	def writeflat(self,flatdelim='\t'):
		outf = open('jobstats.tsv','w')
		outf.write(flatdelim.join(['JobId']+['time']+sorted(self.fields))+'\n')
		for jobid,job in sorted(self.jobs.items(), key=lambda x: x[0]):
			outf.write('%s\n' %job.flatstr(self.fields,flatdelim))
		outf.close()

	def writejson(self):
		json.dump(self.jobs, open('jobstats.json','w'))

	def plot(self):
		plotfields = [
		#	'Resource_List.mem',
		#	'Resource_List.ncpus',
			'resources_used.cpupercent',
		#	'resources_used.cput',
			'resources_used.mem',
		#	'resources_used.ncpus',
		#	'resources_used.vmem',
		]

		fig,axes = pyplot.subplots(len(plotfields),1)
		#fig.set_size_inches(6,6)
		rmarg = 0.7
		pyplot.subplots_adjust(hspace=0.3,bottom=0.1,top=0.95,right=rmarg)

		jobids = sorted(self.jobs.keys())

		cmap = matplotlib.cm.get_cmap('tab20')
		ncolors = 20
		def coli(jobid):
			return(float(jobids.index(jobid) % ncolors)/ncolors)

		# currently system time is in seconds, convert to minutes
		def tfmt(t):
			return(float(int(t)-int(self.tstart))/60)

		for i in range(len(plotfields)):
			fld = plotfields[i]
			maxy = 0
			for jobid,job in sorted(self.jobs.items(), key=lambda x: x[0]):
				x = []
				y = []
				for t,stats in sorted(job.items(), key=lambda x: int(x[0])):
					x.append(tfmt(t))
					y.append(float(stats[fld]))
				axes[i].plot(x,y,label=jobid,color=cmap(coli(jobid)))
				maxy = max(maxy,max(y))

				# add a dashed line for the requested/allocated memory for this job
				if fld in ['resources_used.mem','resources_used.vmem']:
					x = []
					y = []
					for t,stats in sorted(job.items(), key=lambda x: int(x[0])):
						x.append(tfmt(t))
						y.append(float(stats['Resource_List.mem']))
					axes[i].plot(x,y,label=jobid,color=cmap(coli(jobid)),dashes=(2,2))
					maxy = max(maxy,max(y))

				# add a dashed line for the requested/allocated ncpus for this job
				if fld in ['resources_used.ncpus']:
					x = []
					y = []
					for t,stats in sorted(job.items(), key=lambda x: int(x[0])):
						x.append(tfmt(t))
						y.append(float(stats['Resource_List.ncpus']))
					axes[i].plot(x,y,label=jobid,color=cmap(coli(jobid)),dashes=(2,2))
					maxy = max(maxy,max(y))

				# add a dashed line for the max 'cpupercent' available for the requested/allocated ncpus
				if fld in ['resources_used.cpupercent']:
					x = []
					y = []
					for t,stats in sorted(job.items(), key=lambda x: int(x[0])):
						x.append(tfmt(t))
						y.append(float(stats['Resource_List.ncpus'])*100)
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

if __name__ == "__main__":
	logger = JobLogger()
	logger.log(10,60)
