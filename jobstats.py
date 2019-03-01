#!/usr/bin/env python3

# parses PBSPro qstat -f info for simplied and/or flat job stats

import sys
import time

class Job(dict):
	def __init__(self,jobid):
		dict.__init__(self)
		self.jobid = jobid

	def __str__(self,delim=' = '):
		ll = ['Job Id: %s' %self.jobid]
		for k,v in sorted(self.items(), key=lambda x: x[0]):
			if k == 'jobid': continue
			ll.append('\t%s%s%s' %(k,delim,v))
		return('\n'.join(ll))

	def flatstr(self,delim='\t'):
		l = [self.jobid]
		for k,v in sorted(self.items(), key=lambda x: x[0]):
			if k == 'jobid': continue
			l.append(v)
		return(delim.join(l))

fields = [
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

flatdelim = '\t'

hdr = open('jobstats.header','w')
hdr.write(flatdelim.join(['time']+['JobId']+sorted(fields)))
hdr.close()

kv_delim = ' = '

jobs = {}
jobid = None
for l in sys.stdin:
	if l.startswith('Job Id'):
		jobid = l.strip().split()[2]
		jobs[jobid] = Job(jobid)
	elif kv_delim in l:
		kv = l.strip().split(kv_delim)
		if kv[0] in fields:
			jobs[jobid][kv[0]] = kv[1]

for jobid,job in sorted(jobs.items(), key=lambda x: x[0]):
#	print(job)
	print(flatdelim.join([str(round(time.time())), job.flatstr(flatdelim)]))
