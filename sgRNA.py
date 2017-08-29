#!/usr/bin/env python

# Justin Ashworth 2017
# Clean/flexible/object-oriented Python sgRNA finder
# picks out sites in query sequence (easy/fast with regex)
# optional Doench et al. 'on-target score' [seprately downloaded]
# uses [separately isntalled!] finders (blastn, casOFFinder) to check for off-target sites in input [genome] sequence
# outputs sites table with annotated matches, GTF files for visual checks (e.g. IGV)

# restriction_sites file format is just multiple lines like "EcoRI GAATTC" (additional file provided)

from shutil import copyfile
import sys,os,re
import subprocess
from subprocess import Popen
from optparse import OptionParser
# total_ordering extends custom __eq__ and __lt__ class operators to be sufficient for sorting
from functools import total_ordering
import pickle

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

# avoid crRNAs with an early RNA polymerase III termination signal (TTTT), as they will not even be transcribed!
re_pol_iii = re.compile('TTTT')

class Nuclease:
	def __init__(self,name,pam,tgtlen,cutpos,high_fidelity_positions,us=4,ds=3):
		self.name = name
		self.pam = pam
		self.tgtlen = tgtlen
		self.cutpos = cutpos
		self.high_fidelity_positions=high_fidelity_positions
		# certain scoring metrics (Doench et al.) include additional sequence
		self.us=us
		self.ds=ds
		self.make_re()
	def make_re(self):
		pam = ntexpand(self.pam)
		# certain scoring metrics (Doench et al.) include additional sequence
		# the lookahead assertion (?=) evidently allows overlapping matches for Python regex (crucial!)
		self.re = re.compile('(?=(([ACGT]{%i})(([ACGT]{%i})%s)([ACGT]{%i})))' %(self.us,self.tgtlen,pam,self.ds))
	def findsites_onestrand(self,searchseq,seqn,strandstr):
		sites = []
		for s in self.re.finditer(searchseq):
			(fullseq,us,seq,cr,ds) = s.groups()
			# paranoid simple look-back to make sure Python regex/lookahead found a real site
			if not re.search(fullseq,searchseq):
				msg('ERROR: regex %s for %s constructed/returned a fictitious full site seq %s [fix regex!]' %(self.re.pattern,self.name,fullseq))
				sys.exit()
			# skip sites with an early RNA polymerase III termination signal (TTTT)
			if re_pol_iii.search(seq):
				msg('matching sequence with early RNA pol III terminator (TTTT) skipped %s %i %i %s' %(seq,s.start(),s.end(),strandstr))
				continue
			site = NucSite(seqn,self.name,seq,cr,s.start(),strandstr)
			# certain scoring metrics (Doench et al.) include additional sequence
			site.us = us
			site.ds = ds
			for r_e in REs:
				if r_e.cuts(cr): site.REsites.append(str(r_e))
			sites.append(site)
		return sites
	def findsites(self,searchseq,seqn):
		sites = []
		# seqn and strand char are passed here merely for internal site naming purposes
		sites.extend( self.findsites_onestrand(searchseq,seqn,'+') )
		sites.extend( self.findsites_onestrand(rvs_comp_str(searchseq),seqn,'-') )
		return sites
	def tolerates(self,match,tgt):
		for p in self.high_fidelity_positions:
			if match[p] != tgt[p]: return False
		return True

### add your nuclease definitions here ###
nucleases = {
	# recent research suggests that the necessary sgRNA is only 18-19bp in length (not 20bp)
	# high_fidelity_positions are the 'seed' sequence positions where mismatches aren't tolerated, counting back from the 3' (PAM) side
	# cutpos is cut site, counting back from 3' (PAM) side
	'Cas9' : Nuclease('Cas9',pam='NGG',tgtlen=20,cutpos=-3,high_fidelity_positions=[-1,-2,-3,-4,-5,-6,-7,-8,-9,-10])
	# add more here as desired
}

@total_ordering
class Site:
	sep = '\t'
	def __init__(self,seq,start,strand):
		self.seq      = seq
		self.start    = start
		self.strand   = strand
	def __len__(self):
		return len(self.seq)
	def __str__(self):
		out = [self.seq,str(self.start),str(self.start+len(self.seq)),self.strand]
		return self.sep.join(out)
	# following operators are for storing & sorting members of this class in lists and dicts
	def __key(self): return (self.seq, self.start, self.strand)
	def __hash__(self): return hash(self.__key())
	def __eq__(self,other): return self.__key() == other.__key()
	def __lt__(self,other):
		if self.start < other.start: return True
		if self.start > other.start: return False
		if self.strand == '-' and other.strand == '+': return True
		if self.strand == '+' and other.strand == '-': return False
		return self.seq < other.seq

@total_ordering
class NucSite(Site):
	def __init__(self,parent,nuclease,seq,sgRNA,start,strand):
		Site.__init__(self,seq,start,strand)
		self.parent   = parent
		self.nuclease = nuclease
		self.sgRNA    = sgRNA
		self.REsites  = []
		self.matches  = []
		self.accepted_matches  = []
	def seqmatch(self,otherseq):
		if self.seq.upper() == otherseq.upper(): return True
		return False
	def name(self):
		return '%s_%s_%i_%s' %(self.parent,self.nuclease,self.start,self.strand)
	def __str__(self):
		out = ['TGT',self.name(),self.parent,self.seq,self.sgRNA,str(self.start),str(self.start+len(self.seq)),self.strand]
		if hasattr(self,'ontarget'):
			out.append('%0.2f' %self.ontarget)
		out.append(','.join(self.REsites))
		return self.sep.join(out)
	# following operators are for storing & sorting members of this class in lists and dicts
	def __key(self): return (self.parent, self.nuclease, self.seq, self.start, self.strand)
	def __hash__(self): return hash(self.__key())
	def __eq__(self,other): return self.__key() == other.__key()
	def __lt__(self,other):
		if self.parent < other.parent: return True
		if self.parent > other.parent: return False
		if self.nuclease < other.nuclease: return True
		if self.nuclease > other.nuclease: return False
		if self.start < other.start: return True
		if self.start > other.start: return False
		if self.strand == '-' and other.strand == '+': return True
		if self.strand == '+' and other.strand == '-': return False
		return self.seq < other.seq

# this line instructs compiler to autoexpand __lt__ and __eq__ into full set of comparison operators
@total_ordering
class SiteMatch(Site):
	def __init__(self,gnm,seq,start,strand):
		Site.__init__(self,seq,start,strand)
		self.gnm       = gnm
		self.mismatches = -1
		self.tolerated = False
		self.intronic  = False
	def __str__(self):
		out = [self.gnm,self.seq,str(self.start),str(self.start+len(self.seq)),self.strand,str(self.mismatches)]
		if self.tolerated: out.append('***TOLERATED***')
		else: out.append('')
		if self.intronic: out.append('***INTRONIC***')
		else: out.append('')
		return '\t'.join(out)
	# following operators are for storing & sorting members of this class in lists and dicts
	def __key(self): return (self.gnm, self.seq, self.start, self.strand, self.mismatches, self.tolerated, self.intronic)
	def __hash__(self): return hash(self.__key())
	def __eq__(self,other): return self.__key() == other.__key()
	def __lt__(self,other):
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

	def get_matches(self,nucsites,nucleases,maxmis,flags):
		for nuc in nucleases:
			inputfilename = 'casOFFinder.input.%s' %nuc
			inputfile = open(inputfilename,'w')
			# dir containing fasta files
			inputfile.write('%s\n' %self.fastadir)
			tgt = ''.join(['N' for i in range(nucleases[nuc].tgtlen)]) + nucleases[nuc].pam
			inputfile.write('%s\n' %tgt)
			for site in nucsites:
				inputfile.write('%s %i\n' %(nucsites[site].seq,maxmis))
			inputfile.close()

			outfilename = '%s.results' %inputfilename
			cmd = '%s %s C %s' %(self.exe,inputfilename,outfilename)
			msg(cmd)
			prc = Popen(cmd,stdout=subprocess.PIPE,shell=True)
			prc.wait()

			msg('parsing cas-offinder matches')
			for l in open(outfilename):
				f = l.strip().split()
				seq = f[0]
				gnm = f[1]
				gstart  = int(f[2])
				gseq = f[3]
				gstrand = f[4]

				# cas-offinder doesn't keep track of site names so have to re-match
				for sitename in nucsites:
					if not nucsites[sitename].seqmatch(seq): continue
					nucsites[sitename].matches.append( SiteMatch(gnm,gseq,gstart,gstrand) )
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

	def get_matches(self,nucsites,nuclease,maxmis,flags):
		seqs = []
		for site in nucsites:
			seqs.append( (nucsites[site].name(), nucsites[site].seq ))
		sitesf = 'siteseqs.fa'
		writefasta(seqs,sitesf)

		cmd = '%s %s %s -db %s -query %s' %(self.exe,' '.join(self.flags),' '.join(flags),self.fastafiletosearch,sitesf)
		msg(cmd)
		prc = Popen(cmd,stdout=subprocess.PIPE,shell=True)
		for h in prc.stdout:
			f = h.strip().split()
			sitename = f[0]
			gnm = f[1]
			gseq = f[2]
			site = nucsites[sitename]
			gstart  = int(f[3])
			gend    = int(f[4])
			gstrand = f[5]
			if gstrand=='plus': gstrand='+'
			if gstrand=='minus':
				gstrand='-'
				# blast flips start and end for minus matches: fix so start is always 5' on forward strand genome coordinate
				gstart = gend
			nucsites[sitename].matches.append( SiteMatch(gnm,gseq,gstart,gstrand) )
		msg('done parse')

class Scorer:
	def __init__(self,exe,flags):
		self.exe = exe
		self.flags = flags

class CFDScorer(Scorer):
	testsite = ''.join(['A' for i in range(0,24)]) + 'GGGAAA'
	def __init__(self,exe,flags=[]):
		Scorer.__init__(self,exe,flags)
	def exists(self):
		test = '%s --seq %s' %(exe,self.testsite)
		sys.stderr.write('test %s' %test)
		p=Popen(test,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
		p.wait()
		if p.returncode==0:
			score = float(p.stdout.readline().strip().split()[4])
			sys.stderr.write(' %f (pass)\n' %score)
			return True
		else:
			msg('cfd scorer script %s not found or failed' %exe)
			msg(p.stdout.read())
			msg(p.stderr.read())
			return False
	def score(self,site):
		cmd = '%s --seq %s' %(self.exe,''.join([site.us,site.seq,site.ds]))
		sys.stderr.write(cmd)
		p = Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
		p.wait()
		if p.returncode != 0:
			msg(p.stdout.read())
			msg(p.stderr.read())
		score = float(p.stdout.readline().strip().split()[4])
		sys.stderr.write(' %f\n' %score)
		return score

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
	op.add_option('-o','--ontarget',default='',help='algorithm for on-target score (e.g. CFD if you\'ve downloaded it)')
	op.add_option('-e','--evalue',type=float,default=1000,help='blastn evalue (approx: 1=0-2bp, 100=0-3bp, 1000=0-4bp mismatches...)')
	op.add_option('--maxmis',type=int,default=4,help='maximum mismatches to include when finding off-target sites')
	opt,args = op.parse_args()

	if(len(args)==0):
		print('provide a target sequence file in fasta file format, and a genome/background file in fasta format:\n./sgRNA.py -a blast -r restriction_enzymes -g genome.fa TARGET.fa\n')
		op.print_help()
		sys.exit()

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

	seqs = loadfastas(args)

	# fill name-lookup-able dict of nuc sites, as random lookup will be needed sometimes (e.g. parsing site-finder output files)
	nucsites = {}
	for seqn in sorted(seqs):
		for nuc in nucleases.values():
			for s in nuc.findsites(seqs[seqn],seqn):
				msg(str(s))
				nucsites[s.name()] = s

	pickle.dump(nucsites,open('nucsites.p','wb'))

	if opt.ontarget == 'cfd':
		exe = 'rs2_score_calculator.py'
		scorer = CFDScorer(exe)
		if scorer.exists():
			for s in nucsites:
				if not len(nucsites[s].seq)==23 or nucsites[s].nuclease != 'Cas9':
					msg('CFDScorer: skipping site %s because CFD score is trained only for 20-mer crRNA-PAM Cas9 sequences [ACGT]{4}[ACGT]{20}NGG[ACGT]{3}')
					continue
				nucsites[s].ontarget = scorer.score(nucsites[s])

	# check for off-target sites using finder of choice
	finder = finders[opt.algorithm]
	finder.initialize([opt.genomefile])
	flags = ['-evalue %f' %opt.evalue]
	finder.get_matches(nucsites,nucleases,opt.maxmis,flags)

	# check/flag matches for nuclease tolerance
	# check/flag for cuts in intronic sequence
	# count mismatches, sort
	for s in nucsites:
		site = nucsites[s]
		nuclease = nucleases[site.nuclease]
		ntol = 0
		nclose = 0
		nucsites[s].unique=False
		for i in range(len(site.matches)):
			m = site.matches[i]
			m.mismatches = mismatch_count(m.seq,site.seq)
			if m.mismatches > 0 and m.mismatches < 3: nclose+=1
			if nuclease.tolerates(m.seq, site.seq):
				# this has to index fully into original dict due to the way that Python doesn't use pointers/references
				nucsites[s].matches[i].tolerated=True
				ntol+=1

			# flag for intronic cut site
			cut = m.start + len(m.seq) + nuclease.cutpos
			if m.strand == '-': cut = m.start+3
			if m.gnm in introns:
				for start,end in introns[m.gnm]:
					if cut >= start and cut <= end:
						nucsites[s].matches[i].intronic=True
		if ntol==1 and nclose==0: nucsites[s].unique=True

	outroot = 'sites.%s.%s.%i' %(opt.algorithm,opt.genomefile,opt.maxmis)
	gtf = open('%s.gtf'%outroot,'w')
	gtfu = open('%s.unique.gtf'%outroot,'w')
	for s in sorted(nucsites):
		site = nucsites[s]
		print(str(site))
		for m in sorted(site.matches):
			print('\t'.join(['GNM',site.name(),str(m)]))
			# GTF output for genome-verified sites
			if m.mismatches>0: continue
			gtfline = [m.gnm,'JA',site.nuclease,str(m.start),str(m.start+len(m.seq)),'.',m.strand,'.','parent %s; id %s; seq %s; sgRNA %s; REsites %s' %(site.name(),site.name(),site.seq,site.sgRNA,','.join(site.REsites))]
			if hasattr(m,'ontarget'): gtfline.append('; ontarget %s' %site.ontarget)
			gtf.write('%s\n' %'\t'.join(gtfline))
			if site.unique: gtfu.write('%s\n' %gtfline)
	gtf.close()
	gtfu.close()
	msg('done')
