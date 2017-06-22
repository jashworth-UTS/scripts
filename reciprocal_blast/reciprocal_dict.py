#!/usr/bin/env python

class Model:
	def __init__(self,org,id):
		self.org = org
		self.id = id
	def __str__(self):
		return('%s: %s' %(self.org,self.id))
	def __repr__(self):
		return(self.__str__())

def homologs(homdict,x):
	res = []
	if not type(x) is list: hms = [x]
	else: hms = x
	for hm in hms:
		org = hm.org
		if not org in homdict: print('org %s not present')
		res.append([i for i in homdict[org][hm.id].values()])
	if len(res) > 1: return(res)
	return(res[0])

def best_oneway_dict(fs):
	orgs = {}
	for f in fs:
		#	"Fracy.Arath.best"
		s = f.strip().split('.')
		a = s[0]
		b = s[1]
		print('%s %s' %(a,b))
		if not a in orgs: orgs[a] = {}
		for l in open(f):
			s = l.strip().split()
			ga = s[0]
			gb = s[1]
			if not ga in orgs[a]: orgs[a][ga] = {}
			orgs[a][ga][b] = Model(b,gb)

	print('orgs: %s' %(' '.join(orgs.keys())))
	return(orgs)


def make_dict(fs):
	orgs = {}
	for f in fs:
		s = f.strip().split('.')
		a = s[0]
		b = s[1]
		print('%s %s' %(a,b))
		if not a in orgs: orgs[a] = {}
		if not b in orgs: orgs[b] = {}
		for l in open(f):
			s = l.strip().split()
			ga = s[0]
			gb = s[1]
			if not ga in orgs[a]: orgs[a][ga] = {}
			orgs[a][ga][b] = Model(b,gb)
			if not gb in orgs[b]: orgs[b][gb] = {}
			orgs[b][gb][a] = Model(a,ga)

	print('orgs: %s' %(' '.join(orgs.keys())))
	return(orgs)

#def make_tables(orgs):
	# to do: make org-centric tables in pandas format for usage and writing to file
	# or add to sqlite3
