#!/usr/bin/env python

# Justin Ashworth 2017
# Clean/flexible/object-oriented Python sgRNA finder
# picks out sites in query sequence (easy/fast with regex)
# optional Doench et al. 'on-target score' [seprately downloaded]
# uses [separately isntalled!] finders (blastn, casOFFinder) to check for off-target sites in input [genome] sequence
# outputs sites table with annotated matches, GTF files for visual checks (e.g. IGV)

# restriction_sites file format is just multiple lines like "EcoRI GAATTC" (additional file provided)

import sys,os,re
# total_ordering extends custom __eq__ and __lt__ class operators
from functools import total_ordering
import pickle
import subprocess as sp
import multiprocessing as mp

def msg(_msg): sys.stderr.write('%s\n' %_msg)

comp = {
	'A':'T', 'C':'G', 'G':'C', 'T':'A',
  'R':'Y', 'Y':'R', 'S':'W', 'W':'S',
  'K':'M', 'M':'K', 'B':'V', 'D':'H',
  'H':'D', 'V':'B', 'N':'N', '-':'-', 'U':'A'
}

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
	def __init__(self,name,pam,crlen,pampos,cutpos,us=0,ds=0):
		self.name = name
		self.pam = pam
		self.pampos = pampos
		self.crlen = crlen
		self.cutpos = cutpos
		# certain scoring metrics (Doench et al.) include additional sequence
		self.us=us
		self.ds=ds
		self.make_re()
	def pattern(self):
		site = ['N' for i in range(len(self.pam)+self.crlen)]
		for i in range(len(self.pam)): site[self.pampos+i] = self.pam[i]
		return ''.join(site)
	def make_re(self):
		site = self.pattern()
		# the lookahead assertion (?=) evidently allows overlapping matches for Python regex (crucial!)
		site = ntexpand(site)
		# certain scoring metrics (Doench et al.) include additional sequence
		self.re = re.compile('(?=([ACGT]{%i})(%s)([ACGT]{%i}))' %(self.us,site,self.ds))
	def findsites_onestrand(self,searchseq,seqn,strandstr):
		sites = []
		seql = len(searchseq)
		for s in self.re.finditer(searchseq):
			(us,siteseq,ds) = s.groups()
			# skip sites with an early RNA polymerase III termination signal (TTTT)
			if re_pol_iii.search(siteseq):
				msg('matching sequence with early RNA pol III terminator (TTTT) skipped %s %i %i %s' %(siteseq,s.start(),s.end(),strandstr))
				continue
			# crRNA sequence defined as all site letters not in the PAM
			cr = []
			for i in range(len(siteseq)):
				if i >= self.pampos and i < self.pampos + len(self.pam): continue
				cr.append(siteseq[i])
			cr = ''.join(cr)
			start = s.start()
			if strandstr == '-':
				start = seql - start
			site = NucSite(seqn,self.name,us,siteseq,cr,ds,start,strandstr)
			# certain scoring metrics (Doench et al.) include additional sequence
			for r_e in REs:
				if r_e.cuts(siteseq): site.REsites.append(str(r_e))
			sites.append(site)
		return sites
	def findsites(self,searchseq,seqn):
		sites = []
		# seqn and strand char are passed here merely for internal site naming purposes
		sites.extend( self.findsites_onestrand(searchseq.upper(),seqn,'+') )
		sites.extend( self.findsites_onestrand(rvs_comp_str(searchseq.upper()),seqn,'-') )
		return sites

### nuclease definitions ###
nucleases = {
	# recent research suggests that the necessary crRNA is only 18-19bp in length (not 20bp)
	# cutpos is cut site, counting back from 3' (PAM) side
	'SpCas9' : Nuclease('SpCas9',pam='NGG',crlen=20,pampos=20,cutpos=-3,us=4,ds=3),
	'LbCpf1' : Nuclease('LbCpf1',pam='TTTV',crlen=20,pampos=0,cutpos=21),
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
	def __init__(self,parent,nuclease,us,seq,crRNA,ds,start,strand):
		Site.__init__(self,seq,start,strand)
		self.parent   = parent
		self.nuclease = nuclease
		self.us       = us
		self.crRNA    = crRNA
		self.ds       = ds
		self.REsites  = []
		self.matches  = []
		self.accepted_matches  = []
	def seqmatch(self,otherseq):
		if self.seq.upper() == otherseq.upper(): return True
		return False
	def crmatch(self,otherseq):
		if self.crRNA.upper() == otherseq.upper(): return True
		return False
	def pamstart(self): return len(self.crRNA)
	def thirtymer(self):
		return '%s%s%s' %(self.us,self.seq,self.ds)
	def name(self):
		return '%s_%s_%i_%s' %(self.parent,self.nuclease,self.start,self.strand)
	def __str__(self):
		out = ['TGT',self.name(),self.parent,self.seq,self.crRNA,str(self.start),str(self.start+len(self.seq)),self.strand]
		if hasattr(self,'ontarget'): out.append('%0.2f' %self.ontarget)
		out.append(','.join(self.REsites))
		return self.sep.join(out)
	# following operators are for storing & sorting members of this class in lists and dicts
	def __key(self): return (self.parent, self.nuclease, self.start, self.strand, self.seq)
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
		self.intronic  = False
	def __str__(self):
		out = [self.gnm,self.seq,str(self.start),str(self.start+len(self.seq)),self.strand,str(self.mismatches)]
		if hasattr(self,'offtarget'): out.append('%0.2f' %self.offtarget)
		if self.intronic: out.append('***INTRONIC***')
		return '\t'.join(out)
	# following operators are for storing & sorting members of this class in lists and dicts
	def __key(self): return (self.gnm, self.seq, self.start, self.strand, self.mismatches, self.intronic)
	def __hash__(self): return hash(self.__key())
	def __eq__(self,other): return self.__key() == other.__key()
	def __lt__(self,other):
		if hasattr(self,'offtarget') and hasattr(other,'offtarget'):
			# higher offtarget scores are most important to consider: descending sort
			if self.offtarget > other.offtarget: return True
			if self.offtarget < other.offtarget: return False
		if self.mismatches < other.mismatches: return True
		if self.mismatches > other.mismatches: return False
		if self.gnm < other.gnm: return True
		if self.gnm > other.gnm: return False
		return self.start < other.start

class CFDScorer:
	def __init__(self,threshold,path='/Users/justin/code/Doench_et_al_CRISPRCas9/CFD_Scoring'):
		self.threshold = threshold
		self.path = path
		self.mm_scores = pickle.load(open('%s/mismatch_score.pkl' %self.path,'rb'))
		self.pam_scores = pickle.load(open('%s/pam_scores.pkl' %self.path,'rb'))
	def works(self):
		return hasattr(self,'mm_scores') and hasattr(self,'pam_scores')

	def score(self,tgt,off):
		ltgt = len(tgt)
		loff = len(off)
		if not ltgt == 23 and ltgt==loff or 'N' in off:
#			msg('WARNING CFDScorer: trained only for 20-mer crRNA-PAM SpCas9 sequences with no N bases (site %s)' %off)
			return 0
		score = 1
		tgt = tgt.replace('T','U')
		off = off.replace('T','U')
		# this number (20) is hardcoded because of the study/scoring matrix authors created specifically for SpCas9
		for i in range(20):
			if off[i] == tgt[i]: continue
			key = 'r'+tgt[i]+':d'+rvs_comp_str(off[i])+','+str(i+1)
			score*= self.mm_scores[key]
		score*= self.pam_scores[off[-2:]]
		return score

class SiteFinder:
	def __init__(self,exe,flags=[]):
		self.exe = exe
		self.flags = flags

class CasOFFinder(SiteFinder):
	def __init__(self,exe,flags):
		SiteFinder.__init__(self,exe,flags)

	def initialize(self,fastafilestosearch,offscorer,flags=[]):
		self.fastafilestosearch = fastafilestosearch
		self.offscorer = offscorer
		self.fastadir = 'fasta'
		if os.path.exists(self.fastadir):
			if not os.path.isdir(self.fastadir):
				msg('\'fasta\' needs to be directory')
				sys.exit()
		else: os.mkdir(self.fastadir)
		from shutil import copyfile
		for f in self.fastafilestosearch:
			copyfile(f,'%s/%s' %(self.fastadir,f))

	def write_inputfile(self,nucsites,nuclease,maxmis,scorethreshold):
		inputfilename = 'casOFFinder.input.%s.%i.%f' %(nuclease.name,maxmis,scorethreshold)
		inputfile = open(inputfilename,'w') # dir containing fasta files
		inputfile.write('%s\n' %self.fastadir)
		# for SpCas9, people are a) searching for NGG sites, but b) using 'NRG' for SpCas9 in casOFFinder these days since this includes known off-targets
		pat = nuclease.pattern()
		if nuclease.name=='SpCas9': pat = pat.replace('NGG','NRG')
		inputfile.write('%s\n' %pat)
		for site in nucsites.values():
			inputfile.write('%s %i %f\n' %(site.seq,maxmis,scorethreshold))
		inputfile.close()
		return inputfilename

	def get_matches(self,nucsites,nuclease,maxmis,scorethreshold):
		infile = self.write_inputfile(nucsites,nuclease,maxmis,scorethreshold)
		if nuclease.name=='SpCas9': self.exe = 'cas-offinder-CFD'
		# passing '-' as outfilename to cas-offinder means STDOUT
		cmd = '%s %s C - pam_scores.tsv cfd_mismatch_scores.tsv' %(self.exe,infile)
		msg(cmd)

		# get proc STDOUT and filter by offtarget score in real time
		proc = sp.Popen(cmd,stdout=sp.PIPE,shell=True)
		noff=0
		ndropped=0
		# implementing threaded CFD scoring/filtering in native cas-offinder c++ would be the best performing solution to this bottleneck for large numbers of sites
		cf=open('CasOFF.output','w')
		for l in proc.stdout:
			cf.write(l)
			noff+=1
			f = l.split()
			seq = f[0]
			cr = seq[:-3]
			gnm = f[1]
			gstart  = int(f[2])
			gseq = f[3]
			gstrand = f[4]
			score = 0
			if len(f)>6: score = float(f[6])
			site = SiteMatch(gnm,gseq,gstart,gstrand)
			site.offtarget=score
			nucsites[seq].matches.append(site)
		cf.close()
		return noff,ndropped

class BlastFinder(SiteFinder):
	def __init__(self,exe,flags):
		SiteFinder.__init__(self,exe,flags)

	def initialize(self,fastafilestosearch,offscorer,flags=[]):
		self.flags.extend(flags)
		self.fastafiletosearch=fastafilestosearch[0]
		if os.path.exists('%s.nsq' %self.fastafiletosearch): return
		makeblastdb = 'makeblastdb -in %s -dbtype nucl -max_file_sz 2GB' %self.fastafiletosearch
		msg(makeblastdb)
		mdb = sp.Popen(makeblastdb,shell=True)
		mdb.wait()

	def get_matches(self,nucsites,nuclease,maxmis):
		seqs = []
		for site in nucsites:
			seqs.append( (nucsites[site].name(), nucsites[site].seq ))
		sitesf = 'siteseqs.fa'
		writefasta(seqs,sitesf)

		msg(self.fastafiletosearch)
		cmd = '%s %s -db %s -query %s' %(self.exe,' '.join(self.flags),self.fastafiletosearch,sitesf)
		msg(cmd)
		prc = sp.Popen(cmd,stdout=sp.PIPE,shell=True)
		noff = 0
		# real-time/concurrent parse of STDOUT(?)
		for h in prc.stdout:
			f = h.strip().split()
			sitename = f[0]
			gnm = f[1]
			qry = f[2]
			gseq = f[3]
			gstart  = int(f[4])
			gend    = int(f[5])
			gstrand = f[6]
			if gstrand=='plus': gstrand='+'
			if gstrand=='minus':
				gstrand='-'
				# blast flips start and end for minus matches: fix so start is always 5' on forward strand genome coordinate
				gstart = gend
			nucsites[qry[:-3]].matches.append( SiteMatch(gnm,gseq,gstart,gstrand) )
			noff+=1
		msg('loaded %i blast matches' %noff)
		return noff,0

class RuleSet2Scorer:
	import model_comparison
	def __init__(self,path='/home/justin/code/Doench_et_al_CRISPRCas9/Rule_Set_2_scoring_v1/saved_models'):
		mf = '%s/%s' %(path,'V3_model_nopos.pickle')
		self.model = pickle.load(open(mf,'rb'))
		repickle = mf + '.repickled'
		pickle.dump(self.model,open(repickle,'wb'))
		msg('RuleSet2Scorer: loaded model from %s.' %mf)
	def works(self): return hasattr(self,'model')
	def score(self,seq):
		if len(seq) != 30:
			msg('RuleSet2Scorer calculates on-target scores for 30nt sequences only.')
			return 0
		if seq[25:27] != 'GG':
			msg('RuleSet2Scorer calculates on-target scores for sgRNAs with NGG PAM only.')
			return 0
		return self.model_comparison.predict(seq, -1, -1, model=self.model)

finders = {
	'blast' : BlastFinder('blastn',[
		'-outfmt "6 qseqid sseqid qseq sseq sstart send sstrand"',
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
	from optparse import OptionParser
	op=OptionParser()
	op.add_option('-g','--genomefile',help='genome/background sequence file in fasta format')
	op.add_option('-r','--refile',default='restriction_enzymes',help='restriction enzymes file in "Name Site" format')
	op.add_option('-m','--models',default='',help='gene models file in GFF/GTF format (must define intron elements')
	op.add_option('-a','--algorithm',default='blast',help='algorithm for finding off-target sites; options: blast, casoff')
	op.add_option('-o','--ontarget',default='',help='algorithm for on-target score (e.g. Doench Rule Set 2 if you\'ve downloaded/configed it)')
	op.add_option('-f','--offtarget',default='',help='algorithm for off-target score (e.g. Doench CFD if you\'ve downloaded/configed it)')
	op.add_option('-t','--offtarget_threshold',type='float',default=0.2,help='if an off-target has a score and it\'s below this number, don\'t include it in output lists')
	op.add_option('-e','--evalue',type=float,default=1000,help='blastn evalue (approx: 1=0-2bp, 100=0-3bp, 1000=0-4bp mismatches...)')
	op.add_option('--maxmis',type=int,default=4,help='maximum mismatches to include when finding off-target sites')
	opt,args = op.parse_args()

	if(len(args)==0):
		print('provide a target sequence file in fasta file format, and a genome/background file in fasta format:\n./sgRNA.py -a blast -r restriction_enzymes -g genome.fa TARGET.fa\n')
		op.print_help()
		sys.exit()

	# set up the scorers,
	# and multiprocessing pool /before/ loading much of anything into memory--Pool() procs duplicate the memory of the parent process upon creation
	onscorer = None
	if opt.ontarget == 'rs2':
		onscorer = RuleSet2Scorer()
		def onscore(x): return x[0], onscorer.score(x[1])
	offscorer = None
	if opt.offtarget == 'cfd':
		offscorer = CFDScorer(opt.offtarget_threshold)
		# note I've made a custom version of Cas-OFFinder that applies CFD scoring/filtering internally--use that instead if available

	# Pool must be instantiated /after/ the definition of any functions (e.g. above) that Pool::map() will call later
	procpool = mp.Pool()

	nuclease = nucleases['SpCas9']
#	nuclease = nucleases['LbCpf1']
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
	# fill lookup-able dict of nuclease sites: random lookup will be needed sometimes (e.g. parsing cas-offinder output)
	# using site crRNA as dict key
	nucsites = {}
	nonunique = []
	for seqn in sorted(seqs):
		for s in nuclease.findsites(seqs[seqn],seqn):
			msg(str(s))
			if s.crRNA in nucsites:
				# remember non-unique sites for removal later
				nonunique.append(s.crRNA)
				continue
			nucsites[s.seq] = s

	#	pickle.dump(nucsites,open('nucsites.p','wb'))

	for nu in nonunique:
		msg('sites with seq %s omitted because not unique!' %ns)
		nucsites.pop(nu)

	if onscorer and onscorer.works():
		msg('onscores (%s) for %i found nuclease sites' %(onscorer.__class__.__name__,len(nucsites)))
		# multiprocess onscore calculation
		onscores = procpool.map(onscore, [(k,nucsites[k].thirtymer()) for k in nucsites.keys()])
		for key,score in onscores: nucsites[key].ontarget = score

	if offscorer and offscorer.works():
		msg('offscores: %s' %offscorer.__class__.__name__)

	# find off-target sites using finder of choice (probably casOFFinder, but maybe blast or something for comparison purposes)
	finder = finders[opt.algorithm]
	msg(opt.genomefile)
	finder.initialize([opt.genomefile],offscorer,flags=['-evalue %f' %opt.evalue])
	noff,ndropped = finder.get_matches(nucsites,nuclease,opt.maxmis,opt.offtarget_threshold)

	pskipped = 0
	if ndropped > 0:
		pskipped = 100*ndropped/noff
		msg('parsed %i genomic sites, dropped %i (%2.0f%s) with %s score less than %f' %(noff,ndropped,pskipped,'%',offscorer.__class__.__name__,opt.offtarget_threshold))

	msg('re-tallying mismatches')
	for site in nucsites.values():
		for i in range(len(site.matches)):
			m = site.matches[i]
			m.mismatches = mismatch_count(m.seq,site.seq)

	if len(introns) > 0:
		msg('checking for intronic cutsites')
		for site in nucsites.values():
			cutpos = nuclease.cutpos
			for i in range(len(site.matches)):
				# flag for intronic cut site
				cut = m.start + len(m.seq) + cutpos
				if m.strand == '-': cut = m.start+3
				if m.gnm in introns:
					for start,end in introns[m.gnm]:
						if cut >= start and cut <= end:
							site.matches[i].intronic=True

	msg('GTF files')
	outroot = '%s.sites.%s.%s.%i' %(args[0],opt.algorithm,opt.genomefile,opt.maxmis)
	gtf = open('%s.gtf'%outroot,'w')
	for k,site in sorted(nucsites.items(), key=lambda x: x[1]):
		print(str(site))
		for m in sorted(site.matches):
			print('\t'.join(['OFF',site.name(),str(m)]))
			if m.mismatches==0:
				# GTF output of verified on-target sites mapped to [genome] sequence searched
				gtfline = [m.gnm,'JA',site.nuclease,str(m.start),str(m.start+len(m.seq)),'.',m.strand,'.','parent %s; id %s; seq %s; crRNA %s; REsites %s' %(site.name(),site.name(),site.seq,site.crRNA,','.join(site.REsites))]
				if hasattr(site,'ontarget'): gtfline.append('; ontarget %s' %site.ontarget)
				if hasattr(m,'offtarget'): gtfline.append('; offtarget %s' %m.offtarget)
				gtf.write('%s\n' %'\t'.join(gtfline))
			# to do: GTF output for off-target sites w/ genome locations/IGV searchability?
#			else:

	gtf.close()
	msg('done')
