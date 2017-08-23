#!/usr/bin/env python
import sys,re,os
import subprocess
from subprocess import Popen
from optparse import OptionParser

class Nuclease:
	def __init__(self,name,pam,sglength):
		self.name = name
		self.pam = pam
		self.sgl = sglength
		self.make_re()
	def make_re(self):
		self.re = re.compile('(([ACGT]{%i})%s)' %(self.sgl,self.pam))

### add your nuclease definitions here ###
nucleases = {
	'Cas9' : Nuclease('Cas9','[ACGT]GG',20),
#	'another' : Nuclease('name','[ACGT]TT',20),
}

class Site:
	sep = '\t'
	header = sep.join(['parent','nuclease','site','sgRNA','start','end','strand'])
	def __init__(self,parent,nuclease,seq,sgRNA,start,end,strand):
		self.parent   = parent
		self.nuclease = nuclease
		self.seq      = seq
		self.sgRNA    = sgRNA
		self.start    = start
		self.end      = end
		self.strand   = strand
		self.matches  = []
	def __len__(self):
		return len(self.seq)
	def name(self):
		strand = 'fwd'
		if self.strand=='antisense': strand = 'rvs'
		return '%s_%i_%i_%s' %(self.parent,self.start,self.end,strand)
	def __str__(self):
		return self.sep.join([self.name(),self.nuclease,self.seq,self.sgRNA,str(self.start),str(self.end),self.strand])

### util functions ###
def writefasta(seqs,fname):
	f=open(fname,'w')
	for s in seqs:
		if not s: continue
		f.write('>%s\n%s\n' %(s[0],s[1]))
	f.close()

def loadfastas(fastas):
	seqs = {}
	for fasta in fastas:
		if not fasta or not os.path.exists(fasta):
			sys.stderr.write('skipping missing file %s\n' %fasta)
			continue
		seqname=''
		seq=[]
		sys.stderr.write('reading sequence from %s\n' %fasta)
		for line in open(fasta):
			if line.startswith('>'):
				if seqname != '' and len(seq)>0:
					seqs[seqname] = ''.join(seq)
					seq=[]
				seqname=re.sub('>','',line.strip()).strip()
			else:
				seq.append( line.strip() )
		# last hanging seq
		if seqname != '' and len(seq)>0:
			seqs[seqname] = ''.join(seq)
	return seqs

comp = {
	'A':'T', 'C':'G', 'G':'C', 'T':'A',
	'a':'t', 'c':'g', 'g':'c', 't':'a',
  'R':'Y', 'Y':'R', 'S':'W', 'W':'S',
  'K':'M', 'M':'K', 'B':'V', 'D':'H',
  'H':'D', 'V':'B', 'N':'N', }

def compbase(char):
	try: return comp[char]
	except:
		sys.stderr.write('failed to complement base %s\n' %char)
		return char

def rvs_comp_str(seq):
	return(''.join([ compbase(base) for base in reversed(seq) ]))

if __name__ == "__main__":
	op=OptionParser()
	op.add_option('-g','--genomefile')
	opt,args = op.parse_args()

	seqs = loadfastas(args)

	sites = {}
	siteorder = []
	for seqn in sorted(seqs):
		for nuc in nucleases.values():
			for s in nuc.re.finditer(seqs[seqn]):
				site = Site(seqn,nuc.name,s.group(1),s.group(2),s.start(),s.end(),'fwd')
				sites[site.name()] = site
				siteorder.append(site.name())
			for s in nuc.re.finditer(rvs_comp_str(seqs[seqn])):
				site = Site(seqn,nuc.name,s.group(1),s.group(2),s.start(),s.end(),'rvs')
				sites[site.name()] = site
				siteorder.append(site.name())

#	print(Site.header)
	seqs = []
	for site in siteorder:
#		print(str(site))
		seqs.append( (sites[site].name(), sites[site].seq ))
	sitesf = 'siteseqs.fa'
	writefasta(seqs,sitesf)

	if not os.path.exists('%s.nsq' %opt.genomefile):
		Popen('makeblastdb -in %s -dbtype nucl -max_file_size 2GB')

	# fast genome matching strategy: use blastn to get small kmer-based near-exact matches
	# can re-check later by precise walking method (JA pssm++ program)

	# run blastn
	outfmt = '"6 qseqid sseqid mismatch sseq sstart send sstrand"'
# 	cmd = 'blastn -task blastn-short -query %s -db %s -ungapped -outfmt %s -evalue 1 -qcov_hsp_perc 100 > %s.blastn' %(sitesf,opt.genomefile,outfmt,sitesf)
 	cmd = 'blastn -task blastn-short -query %s -db %s -ungapped -outfmt %s -evalue 1 -qcov_hsp_perc 100' %(sitesf,opt.genomefile,outfmt)
	blastn = Popen(cmd,stdout=subprocess.PIPE,shell=True)

	for h in blastn.stdout:
		f = h.strip().split()
		sitename = f[0]
		gnm = f[1]
		mis = int(f[2])
		gsq = f[3]
		gstart  = int(f[4])
		gend    = int(f[5])
		gstrand = f[6]
		if gstrand=='plus': gstrand='fwd'
		if gstrand=='minus': gstrand='rvs'
		sites[sitename].matches.append( (gnm,mis,gsq,gstart,gend,gstrand) )

	for s in siteorder:
		site = sites[s]
		print(str(site))
		for m in site.matches:
			print('\t'.join([site.name(),m[0],m[2],str(m[1]),str(m[3]),str(m[4]),m[5]]))

