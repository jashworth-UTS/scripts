#!/usr/bin/env python
import sys, os, re, math

bases = ['A','C','G','T','N']

codons = {
	'A' : ['GCA','GCC','GCG','GCT'],
	'C' : ['TGC','TGT'],
	'D' : ['GAC','GAT'],
	'E' : ['GAA','GAG'],
	'F' : ['TTC','TTT'],
	'G' : ['GGA','GGC','GGG','GGT'],
	'H' : ['CAC','CAT'],
	'I' : ['ATA','ATC','ATT'],
	'K' : ['AAA','AAG'],
	'L' : ['CTA','CTC','CTG','CTT','TTA','TTG'],
	'M' : ['ATG'],
	'N' : ['AAC','AAT'],
	'P' : ['CCA','CCC','CCG','CCT'],
	'Q' : ['CAA','CAG',],
	'R' : ['CGA','CGC','CGG','CGT','AGA','AGG'],
	'S' : ['AGC','AGT','TCA','TCC','TCG','TCT'],
	'T' : ['ACA','ACC','ACG','ACT'],
	'V' : ['GTA','GTC','GTG','GTT'],
	'W' : ['TGG'],
	'Y' : ['TAC','TAT'],
	'.' : ['TAA','TAG','TGA'],
}

translation_code = {}
n = 0
for aa,codonlist in codons.items():
	n += len(codonlist)
	for codon in codonlist:
		translation_code[codon] = aa

def msg(_msg): sys.stderr.write('%s\n' %_msg)

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

def all_combinations_gen( pos, length, seq, seqs, choices ):
	if seq == []:
		for i in range(length): seq.append( '' )
	for choice in choices[pos]:
		seq[pos] = choice
		if pos == length-1:
			seqs.append( ''.join( seq ) )
		else: all_combinations_gen( pos+1, length, seq, seqs, choices )

def triplet_count(seqs):
	choices = [ bases for i in range(3) ]
	triplets = []
	all_combinations_gen(0,3,[],triplets,choices)
#	print(triplets)
	counts = {}
	for t in triplets: counts[t] = 0
	for seq in seqs.values():
		i=0
		while i < len(seq)-2:
			# dont bother filtering 'N's here--too slow. Do it later
			counts[seq[i:(i+3)]] += 1
			i += 3
	# now drop any 'N' triplets--fast and easy
	drops = []
	for t in counts:
		if 'N' in t: drops.append(t)
	for d in drops: counts.pop(d)
	return counts

def counts_freq(counts):
	freqs = {}
	total = 0
	for c in counts.values(): total += c
	for t,c in counts.items():
		freqs[t] = float(c)/total
	return freqs

def counts_summary(counts):
	freqs = counts_freq(counts)
	for t,fq in freqs.items():
		print('%s %0.4f' %(t,fq))

def compare(qrf,bgf):
	msg('AA-specific codon frequencies comparison:')
	qry = loadfastas(qrf)
	bg = loadfastas(bgf)
	counts = triplet_count(bg)
	for n,q in qry.items():
		i = 0
		sumllr = 0
		while i < len(q)-2:
			cod = q[i:(i+3)]
			aa = translation_code[cod]
			aatot = 0
			opt = ''
			optcnt = 0
			for allcod in codons[aa]:
				cnt = counts[allcod]
				aatot += cnt
				if cnt > optcnt:
					opt = allcod
					optcnt = cnt
			fq = float(counts[cod])/aatot
			fqopt = float(optcnt)/aatot
			lr = fq/fqopt
			flag = '    '
			if lr < 1:   flag = '*   '
			if lr < 0.7: flag = '**  '
			if lr < 0.5: flag = '*** '
			if lr < 0.3: flag = '****'
			sumllr += math.log(lr)
			print('%i\t%s\t%0.3f\t%s\t%s\t%0.3f\t%s\t%0.2f' %(i,cod,fq,aa,opt,fqopt,flag,lr))
			i+=3
		print('sum log liklihood ratio: %g' %sumllr)

def codefreq(fnames):
	seqs = loadfastas(fnames)
	counts = triplet_count(seqs)
	counts_summary(counts)

def do_codefreq():
	codefreq(sys.argv[1:])

def do_compare():
	qrf = sys.argv[1]
	bg = sys.argv[2]
	compare([qrf],[bg])

if __name__  == "__main__":
	do_compare()
