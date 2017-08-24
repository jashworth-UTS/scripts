#!/usr/bin/env python

# JA 2017
# simple sgRNA finder for fasta input target sequences
# uses blastn to check for off-target sites in input [genome] sequence
# outputs sites table with blastn matches, and GTF files for visual checks (e.g. IGV)

# restriction_sites file format is just multiple lines like "EcoRI GAATTC"
# (additional file is provided for this)

import sys,re,os
import subprocess
from subprocess import Popen
from optparse import OptionParser

def msg(_msg): sys.stderr.write('%s\n' %_msg)

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

comp = {
	'A':'T', 'C':'G', 'G':'C', 'T':'A',
  'R':'Y', 'Y':'R', 'S':'W', 'W':'S',
  'K':'M', 'M':'K', 'B':'V', 'D':'H',
  'H':'D', 'V':'B', 'N':'N', '-':'-' }

def compbase(char):
	try: return comp[char]
	except:
		msg('failed to complement base %s' %char)
		return char

def rvs_comp_str(seq):
	return(''.join([ compbase(base) for base in reversed(seq) ]))

class REnz:
	def __init__(self,name='',site=''):
		self.name = name
		self.site = site
		self.re = re.compile(site)
	def cuts(self,seq):
		if self.re.search(seq): return True
		if self.re.search(rvs_comp_str(seq)): return True
		return False
	def __str__(self):
		return '%s(%s)' %(self.name,self.site)

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
		self.REsites  = []
		self.matches  = []
		self.gnm      = None
	def __len__(self):
		return len(self.seq)
	def name(self):
		strand = 'fwd'
		if self.strand=='antisense': strand = 'rvs'
		return '%s_%i_%i_%s' %(self.parent,self.start,self.end,strand)
	def __str__(self):
		out = ['TGT',self.name(),self.nuclease,self.seq,self.sgRNA,str(self.start),str(self.end),self.strand]
		out.append(','.join(self.REsites))
		return self.sep.join(out)

### util functions ###
def read_restriction_sites(fname):
	msg('reading %s' %fname)
	res = []
	for l in open(fname):
		fs = l.strip().split()
		res.append(REnz(fs[0],fs[1]))
	return res

def writefasta(seqs,fname):
	msg('writing %s' %fname)
	f=open(fname,'w')
	for s in seqs:
		if not s: continue
		f.write('>%s\n%s\n' %(s[0],s[1]))
	f.close()

def loadfastas(fastas):
	seqs = {}
	for fasta in fastas:
		if not fasta or not os.path.exists(fasta):
			msg('skipping missing file %s' %fasta)
			continue
		msg('reading %s' %fasta)
		seqname=''
		seq=[]
		msg('reading sequence from %s' %fasta)
		for line in open(fasta):
			if line.startswith('>'):
				if seqname != '' and len(seq)>0:
					seqs[seqname] = ''.join(seq).upper()
					seq=[]
				seqname=re.sub('>','',line.strip()).strip()
			else:
				seq.append( line.strip() )
		# last hanging seq
		if seqname != '' and len(seq)>0:
			seqs[seqname] = ''.join(seq).upper()
	return seqs

if __name__ == "__main__":
	op=OptionParser()
	op.add_option('-g','--genomefile',help='genome/background sequence file in fasta format')
	op.add_option('-r','--refile',default='restriction_enzymes',help='restriction enzymes file in "Name Site" format')
	op.add_option('-m','--models',default='',help='gene models file in GFF/GTF format (must define intron elements')
	op.add_option('-e','--evalue',type=float,default=100,help='blastn evalue (approx: 1=0-2bp, 100=0-3bp, 1000=0-4bp mismatches...)')
	op.add_option('-c','--coverage',type=float,default=90,help='blastn query coverage for off-target search (100%==full-length sites)')
	op.add_option('-u',dest='ungapped',default=False,action='store_true',help='force blastn ungapped alignment')
	opt,args = op.parse_args()

	REs = read_restriction_sites(opt.refile)

	introns = {}
	re_intron = re.compile("(\w+)\s\w+\sintron\s(\d+)\s(\d+)")
	if not opt.models=='':
		msg('reading introns from %s' %opt.models)
		for match in re_intron.finditer(open(opt.models).read()):
			(seqn,start,end) = match.groups()
			if not seqn in introns: introns[seqn] = []
			introns[seqn].append( (int(start),int(end)) )
		msg('read introns for %i contigs' %len(introns))

	if(len(args)==0):
		print('''
provide a target sequence file in fasta file format, and a genome/background file in fasta format:\n
./sgRNA.py -r restriction_enzymes -g genome.fa TARGET.fa\n
%s
		       ''' %opt.help())
		sys.exit()

	seqs = loadfastas(args)

	sites = {}
	siteorder = []
	for seqn in sorted(seqs):
		for nuc in nucleases.values():
			for s in nuc.re.finditer(seqs[seqn]):
				site = Site(seqn,nuc.name,s.group(1),s.group(2),s.start(),s.end(),'fwd')
				for r_e in REs:
					if r_e.cuts(s.group(1)): site.REsites.append(str(r_e))
				sites[site.name()] = site
				siteorder.append(site.name())
			for s in nuc.re.finditer(rvs_comp_str(seqs[seqn])):
				site = Site(seqn,nuc.name,s.group(1),s.group(2),s.start(),s.end(),'rvs')
				for r_e in REs:
					if r_e.cuts(s.group(1)): site.REsites.append(str(r_e))
				sites[site.name()] = site
				siteorder.append(site.name())

	seqs = []
	for site in siteorder:
		seqs.append( (sites[site].name(), sites[site].seq ))
	sitesf = 'siteseqs.fa'
	writefasta(seqs,sitesf)

	if not os.path.exists('%s.nsq' %opt.genomefile):
		makeblastdb = 'makeblastdb -in %s -dbtype nucl -max_file_sz 2GB' %opt.genomefile
		msg(makeblastdb)
		mdb = Popen(makeblastdb,shell=True)
		mdb.wait()
		msg('done makeblastdb')

	# fast genome matching strategy: use blastn to get small kmer-based near-exact matches
	# can re-check later by precise walking method (JA pssm++ program)

	# run blastn
	outfmt = '"6 qseqid sseqid mismatch sseq sstart send sstrand"'
# 	cmd = 'blastn -task blastn-short -query %s -db %s -ungapped -outfmt %s -evalue 1 -qcov_hsp_perc 100 > %s.blastn' %(sitesf,opt.genomefile,outfmt,sitesf)
	ungapped = ''
	if opt.ungapped: ungapped = '-ungapped'
 	cmd = 'blastn -evalue %f -word_size 7 -qcov_hsp_perc %f %s -query %s -db %s -outfmt %s' %(opt.evalue,opt.coverage,ungapped,sitesf,opt.genomefile,outfmt)
	msg(cmd)
	blastn = Popen(cmd,stdout=subprocess.PIPE,shell=True)

	msg('parse blastn...')
	for h in blastn.stdout:
		f = h.strip().split()
		sitename = f[0]
		gnm = f[1]
		mis = int(f[2])
		gsq = f[3]
		gstart  = int(f[4])
		gend    = int(f[5])
		gstrand = f[6]
		cut = gend-3
		if gstrand=='plus': gstrand='+'
		if gstrand=='minus':
			gstrand='-'
			gsq = rvs_comp_str(gsq)
			temp = gstart
			# blast flips start end end for minus matches
			gstart = gend
			gend = temp
			cut = gstart+3
		intron = ''
		if gnm in introns:
			for start,end in introns[gnm]:
				if cut >= start and cut <= end: intron='******INTRON******'
		if not sites[sitename].gnm: sites[sitename].gnm = gnm
		sites[sitename].matches.append( (gnm,gsq,mis,gstart,gend,gstrand,intron) )
	msg('done parse')

	gtff = 'sites.gtf'
	gtf = open(gtff,'w')
	gtfuf = 'unique_sites.gtf'
	gtfu = open(gtfuf,'w')
	for s in siteorder:
		site = sites[s]
		print(str(site))
		for m in site.matches:
			print('\t'.join(['GNM',site.name(),m[0],m[1],str(m[2]),str(m[3]),str(m[4]),m[5],m[6]]))
			# output to GTF only if site is unique
			gtf.write('%s\n' %'\t'.join([m[0],'JA',site.nuclease,str(m[3]),str(m[4]),'.',m[5],'.','parent %s; id %s; seq %s; sgRNA %s; REsites %s' %(site.name(),site.name(),site.seq,site.sgRNA,','.join(site.REsites))]))
			if len(site.matches)==1:
				gtfu.write('%s\n' %'\t'.join([m[0],'JA',site.nuclease,str(m[3]),str(m[4]),'.',m[5],'.','parent %s; id %s; seq %s; sgRNA %s; REsites %s' %(site.name(),site.name(),site.seq,site.sgRNA,','.join(site.REsites))]))
	gtf.close()
	gtfu.close()
	msg('done')
