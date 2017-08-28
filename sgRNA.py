#!/usr/bin/env python

# Justin Ashworth 2017
# Clean/flexible/object-oriented Python sgRNA finder
# uses [serparated isntalled!] finders (blastn, casOFFinder) to check for off-target sites in input [genome] sequence
# outputs sites table with annotated matches, GTF files for visual checks (e.g. IGV)

# restriction_sites file format is just multiple lines like "EcoRI GAATTC" (additional file provided)

from shutil import copyfile
import sys,re,os
import subprocess
from subprocess import Popen
from optparse import OptionParser
from functools import total_ordering

def msg(_msg): sys.stderr.write('%s\n' %_msg)

comp = {
	'A':'T', 'C':'G', 'G':'C', 'T':'A',
  'R':'Y', 'Y':'R', 'S':'W', 'W':'S',
  'K':'M', 'M':'K', 'B':'V', 'D':'H',
  'H':'D', 'V':'B', 'N':'N', '-':'-' }

degen_nucs = {
	'A' : ['A'],
	'C' : ['C'],
	'G' : ['G'],
	'T' : ['T'],
	'R' : ['A','G'],
	'Y' : ['C','T'],
	'M' : ['A','C'],
	'K' : ['G','T'],
	'S' : ['C','G'],
	'W' : ['A','T'],
	'B' : ['C','G','T'],
	'D' : ['A','G','T'],
	'H' : ['A','C','T'],
	'V' : ['A','C','G'],
	'N' : ['A','C','G','T'],
}

def compbase(char):
	try: return comp[char]
	except:
		msg('failed to complement base %s' %char)
		return char

def ntexpand(nts):
	ret = []
	for nt in nts:
		deg = degen_nucs[nt]
		ret.append('[%s]' %''.join(deg))
	return ''.join(ret)

def mismatch_count(a,b):
	cnt = 0
	for i in range(len(a)):
		if a[i] != b[i]: cnt+=1
	return cnt

class Nuclease:
	def __init__(self,name,pam,tgtlen,cutpos,high_fidelity_positions=[]):
		self.name = name
		self.pam = pam
		self.tgtlen = tgtlen
		self.cutpos = cutpos
		self.high_fidelity_positions=high_fidelity_positions
		self.make_re()
	def make_re(self):
		pam = ntexpand(self.pam)
		self.re = re.compile('(([ACGT]{%i})%s)' %(self.tgtlen,pam))
	def tolerates(self,match,tgt):
		for p in self.high_fidelity_positions:
			if match[p] != tgt[p]: return False
		return True

### add your nuclease definitions here ###
nucleases = {
	# recent research suggests that the necessary sgRNA is only 18-19bp in length (not 20bp)
	# high_fidelity_positions are the 'seed' sequence positions where mismatches aren't tolerated, counting back from the 3' (PAM) side
	# cutpos is cut site, counting back from 3' (PAM) side
	'Cas9' : Nuclease('Cas9',pam='NGG',tgtlen=18,cutpos=-3,high_fidelity_positions=[-1,-2,-3,-4,-5,-6,-7,-8,-9,-10])
	# add more here as desired
}

@total_ordering
class SiteMatch:
	def __init__(self,gnm,sq,start,end,strand):
		self.gnm       = gnm
		self.sq        = sq
		self.start     = start
		self.end       = end
		self.strand    = strand
		self.mismatches = -1
		self.tolerated = False
		self.intronic  = False
	def __str__(self):
		out = [self.gnm,self.sq,str(self.start),str(self.end),self.strand,str(self.mismatches)]

		if self.tolerated: out.append('tolerated')
		else: out.append('untolerated')

		if self.intronic: out.append('***INTRONIC***')
		else: out.append('')

		return '\t'.join(out)

	def __eq__(self,other):
		return self.gnm == other.gnm and self.sq == other.sq and self.start == other.start and self.strand == other.strand

	def __lt__(self,other):
	# note that in Python3 they've decided to try removing the __cmp__ way of sorting classes, with no obvious way of supporting the following kind of multi-conditional sorting logic
		if self.tolerated and not other.tolerated: return True
		if not self.tolerated and other.tolerated: return False
		if self.mismatches < other.mismatches: return True
		if self.mismatches > other.mismatches: return False
		if self.gnm < other.gnm: return True
		if self.gnm > other.gnm: return False
		return self.start < other.start

class SiteFinder:
	def __init__(self,exe,flags=[]):
		self.exe = exe
		self.flags = flags

class CasOFFinder(SiteFinder):
	def __init__(self,exe,flags):
		SiteFinder.__init__(self,exe,flags)
	def initialize(self,fastafilestosearch):
		self.fastafilestosearch = fastafilestosearch
		self.fastadir = 'fasta'
		if os.path.exists(self.fastadir):
			if not os.path.isdir(self.fastadir):
				msg('\'fasta\' needs to be directory')
				sys.exit()
		else: os.mkdir(self.fastadir)
		for f in self.fastafilestosearch:
			copyfile(f,'%s/%s' %(self.fastadir,f))

	def get_matches(self,sites,nucleases,maxmis,flags):
		for nuc in nucleases:
			inputfilename = 'casOFFinder.input.%s' %nuc
			inputfile = open(inputfilename,'w')
			# dir containing fasta files
			inputfile.write('%s\n' %self.fastadir)
			tgt = ''.join(['N' for i in range(nucleases[nuc].tgtlen)]) + nucleases[nuc].pam
			inputfile.write('%s\n' %tgt)
			for site in sites:
				inputfile.write('%s %i\n' %(sites[site].seq,maxmis))
			inputfile.close()

			outfilename = '%s.results' %inputfilename
			cmd = '%s %s C %s' %(self.exe,inputfilename,outfilename)
			msg(cmd)
			prc = Popen(cmd,stdout=subprocess.PIPE,shell=True)
			prc.wait()

			for l in open(outfilename):
				f = l.strip().split()
				seq = f[0]
				gnm = f[1]
				gstart  = int(f[2])
				gsq = f[3]
				gstrand = f[4]

				# cas-offinder doesn't keep track of site names so have to re-match
				for sitename in sites:
					if not sites[sitename].seqmatch(seq): continue
					gend = gstart + len(seq)
					sites[sitename].matches.append( SiteMatch(gnm,gsq,gstart,gend,gstrand) )
			msg('done parse')

class BlastFinder(SiteFinder):
	def __init__(self,exe,flags):
		SiteFinder.__init__(self,exe,flags)

	def initialize(self,fastafilestosearch):
		self.fastafiletosearch=fastafilestosearch[0]
		if os.path.exists('%s.nsq' %self.fastafiletosearch): return
		makeblastdb = 'makeblastdb -in %s -dbtype nucl -max_file_sz 2GB' %self.fastafiletosearch
		msg(makeblastdb)
		mdb = Popen(makeblastdb,shell=True)
		mdb.wait()
		msg('done makeblastdb')

	def get_matches(self,sites,nuclease,maxmis,flags):
		seqs = []
		for site in sites:
			seqs.append( (sites[site].name(), sites[site].seq ))
		sitesf = 'siteseqs.fa'
		writefasta(seqs,sitesf)

		cmd = '%s %s %s -db %s -query %s' %(self.exe,' '.join(self.flags),' '.join(flags),self.fastafiletosearch,sitesf)
		msg(cmd)
		prc = Popen(cmd,stdout=subprocess.PIPE,shell=True)
		for h in prc.stdout:
			f = h.strip().split()
			sitename = f[0]
			gnm = f[1]
			gsq = f[2]
			site = sites[sitename]
			gstart  = int(f[3])
			gend    = int(f[4])
			gstrand = f[5]
			if gstrand=='plus': gstrand='+'
			if gstrand=='minus':
				gstrand='-'
#				gsq = rvs_comp_str(gsq)
				temp = gstart
				# blast flips start end end for minus matches
				gstart = gend
				gend = temp
			sites[sitename].matches.append( SiteMatch(gnm,gsq,gstart,gend,gstrand) )
		msg('done parse')

finders = {
	'blast' : BlastFinder('blastn',[
		'-outfmt "6 qseqid sseqid sseq sstart send sstrand"',
		'-word_size 7',
		'-qcov_hsp_perc 100',
		'-ungapped',
	]),
	'casoff' : CasOFFinder('cas-offinder',[]),
}

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
	header = sep.join(['parent','nuclease','site','sgRNA','start','end','strand','mismatches','tolerance','region'])
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
		self.accepted_matches  = []
	def seqmatch(self,seq):
		su = seq.upper()
		if su == self.seq: return True
		if rvs_comp_str(su) == self.seq: return True
		return False
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
	op.add_option('-a','--algorithm',default='blast',help='algorithm for finding off-target sites; options: blast, casoff')
	op.add_option('-e','--evalue',type=float,default=1000,help='blastn evalue (approx: 1=0-2bp, 100=0-3bp, 1000=0-4bp mismatches...)')
	op.add_option('--maxmis',type=int,default=4,help='maximum mismatches to include when finding off-target sites')
	opt,args = op.parse_args()

	if(len(args)==0):
		print('provide a target sequence file in fasta file format, and a genome/background file in fasta format:\n./sgRNA.py -a blast -r restriction_enzymes -g genome.fa TARGET.fa\n')
		op.print_help()
		sys.exit()

	REs = read_restriction_sites(opt.refile)
	finder = finders[opt.algorithm]

	introns = {}
	re_intron = re.compile("(\w+)\s\w+\sintron\s(\d+)\s(\d+)")
	if not opt.models=='':
		msg('reading introns from %s' %opt.models)
		for match in re_intron.finditer(open(opt.models).read()):
			(seqn,start,end) = match.groups()
			if not seqn in introns: introns[seqn] = []
			introns[seqn].append( (int(start),int(end)) )
		msg('read introns for %i contigs' %len(introns))

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

	# check for off-target sites using finder of choice
	finder.initialize([opt.genomefile])
	flags = ['-evalue %f' %opt.evalue]
	finder.get_matches(sites,nucleases,opt.maxmis,flags)

	# check/flag matches for nuclease tolerance
	# check/flag for cuts in intronic sequence
	# count mismatches, sort
	for s in sites:
		site = sites[s]
		nuclease = nucleases[site.nuclease]
		for i in range(len(site.matches)):
			m = site.matches[i]
			m.mismatches = mismatch_count(m.sq,site.seq)
			if nuclease.tolerates(m.sq, site.seq):
				# this has to index fully into original dict due to the way that Python doesn't use pointers/references
				sites[s].matches[i].tolerated=True

			# flag for intronic cut site
			cut = m.end + nuclease.cutpos
			if m.strand == '-': cut = m.start+3
			if m.gnm in introns:
				for start,end in introns[m.gnm]:
					if cut >= start and cut <= end:
						sites[s].matches[i].introniv=True

	gtff = 'sites.gtf'
	gtf = open(gtff,'w')
	gtfuf = 'unique_sites.gtf'
	gtfu = open(gtfuf,'w')
	for s in siteorder:
		site = sites[s]
		print(str(site))
		for m in sorted(site.matches):
			print('\t'.join(['GNM',site.name(),str(m)]))

			# output to GTF only if site is unique
			gtf.write('%s\n' %'\t'.join([m.gnm,'JA',site.nuclease,str(m.start),str(m.end),'.',m.strand,'.','parent %s; id %s; seq %s; sgRNA %s; REsites %s' %(site.name(),site.name(),site.seq,site.sgRNA,','.join(site.REsites))]))
			if len(site.matches)==1:
				gtfu.write('%s\n' %'\t'.join([m.gnm,'JA',site.nuclease,str(m.start),str(m.end),'.',m.strand,'.','parent %s; id %s; seq %s; sgRNA %s; REsites %s' %(site.name(),site.name(),site.seq,site.sgRNA,','.join(site.REsites))]))
	gtf.close()
	gtfu.close()
	msg('done')
