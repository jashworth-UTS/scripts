#!/usr/bin/env python3

import sys,os

# subclass dict to attach methods
class BLAST_nr(dict):
	def __init__(self):
		dict.__init__(self)
	def parseBLAST(self,fname,thresh=100,minalen=50,aggregate=False):
		self.clear()
		# expecting a blast output with -outfmt "6 qseqid qlen sseqid slen pident length evalue bitscore" e.g.:
		# TRINITY_DN197_c5_g1_i1	310	TRINITY_DN197_c5_g1_i1	310	100.000	310	4.86e-163	573

		# collect each query id's longest mate (including self)
		longest_mate = {}
		for l in open(fname):
			ll = l.strip().split()
			pident = float(ll[4])
			if pident < thresh: continue
			alen = int(ll[5])
			if alen < minalen: continue
			qry = ll[0]
			qlen = int(ll[1])
			subj = ll[2]
			slen = int(ll[3])
			if not qry in longest_mate: longest_mate[qry] = (subj,slen)
			if slen > longest_mate[qry][1]: longest_mate[qry] = (subj,slen)
		# longest_mate dict now points to [redundant] set of longest mates (values)

		# transpose to longest mates, with list of 'sub-mates'; keys resultingly non-redudant
		for qry,longest in longest_mate.items():
			subj = longest[0]
			if not subj in self: self[subj] = []
			self[subj].append(qry)

		# unique-ify sub-mates
		for lng in self:
			self[lng] = list(set(self[lng]))

		if aggregate:
			# add exact subclusters into larger metaclusters by referencing sub-mate lists
			# allowing now for inexact matches within combined metaclusters
			topop = []
			for longest,mates in self.items():
				for m in mates:
					# skip self-containment, otherwise will pop itself
					if m==longest: continue
					# for non-self mates--which should all be shorter--add its submates (if existent); later will pop the submate's own list
					if m in self and longest in self:
						currmates = list(self[longest])
						currmates.extend(self[m])
						self[longest] = list(set(currmates))
						topop.append(m)
			sys.stderr.write('popping %i sub-ids after aggregation\n' %len(topop))
			for p in set(topop): self.pop(p)

	def __str__(self):
		out = []
		for lng,mates in self.items():
			out.append('%s:\n\t%s' %(lng,', '.join(mates)))
		return('%s\n' %'\n'.join(out))

	def non_uni_clusters(self):
		nucs = []
		for lng,mates in self.items():
			if len(mates)==1 and lng==mates[0]: continue
			nucs.append( (lng,len(mates),mates) )
		out = []
		for lng,ln,mates in sorted(nucs,key=lambda x: x[1]):
			out.append('%s:\n\t%s' %(lng,', '.join(mates)))
		return('%s\n' %'\n'.join(out))

	# return set of all ids including parent and child ids
	def allids(self):
		ids = []
		for k in self:
			ids.append(k)
			for kk in self[k]: ids.append(kk)
		return(set(ids))

if __name__ == "__main__":
	import sys
	app = BLAST_nr()
	# args are blasttable, minpident, minalen
	bfnm = sys.argv[1]
	minpident = float(sys.argv[2])
	minalnlen = int(sys.argv[3])
	app.parseBLAST(bfnm,minpident,minalnlen)
#	print(app)
#	print(app.non_uni_clusters())
#	print('\n'.join(app.keys()))

	grp_path = 'BLASTn_grps'
	if not os.path.exists(grp_path): os.mkdir(grp_path)
	import Fasta
	fa = Fasta.FastaSeqs()
	fa.loadseqs([sys.argv[4]])
	sfx = '%i.%i' %(int(minpident),minalnlen)
	lngf = open('%s.lng.%s.fa' %(bfnm,sfx),'w')
	for k,mates in app.items():
		lngf.write('>%s\n%s\n' %(k,fa.seqs[k].seq))
		f = open('%s/%s.grp.%s.fa' %(grp_path,k,sfx),'w')
		f.write('>%s\n%s\n' %(k,fa.seqs[k].seq))
		for m in mates:
			f.write('>%s\n%s\n' %(m,fa.seqs[m].seq))
		f.close()
	lngf.close()
	exf = open('%s.extras.%s.fa' %(bfnm,sfx),'w')
	for k in fa.seqs:
		if not k in app:
			exf.write('>%s\n%s\n' %(k,fa.seqs[k].seq))
	exf.close()
