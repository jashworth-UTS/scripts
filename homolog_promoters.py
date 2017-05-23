#!/usr/bin/env python

from diatoms_sqlite import *
import pandas as pd
from AshworthUtil import rvs_comp_str

# is this redundant with diatoms_sqlite import?
con = sq.connect('diatoms.sqlite.db')
c = con.cursor()

# promoter seq region relative to one cds (by cds internal db id)
def promseq(cds,upstream=500,downstream=100):
	(id,org,seq,pid,start,end,strand) = c.execute('select * from cds where id = %i' %cds).fetchone()
	ps = start - upstream
	pe = end + downstream
	if strand == 1:
		ps = start - downstream
		pe = end + upstream

	if ps < 1: ps = 1

	s = c.execute('select seq from genomic where name = "%s" AND org = "%s"' %(seq,org)).fetchone()[0]
	l = len(s)
	if ps >= l:
		print('out of bounds, check feature %s and sequence %s for org %s' %(cds,seq,org))
	if pe > l: pe = l
	s = s[ps-1:pe]
	if strand == 1:
		s = rvs_comp_str(s)
	return(s)

def promoter_for_model(model):
	return ('%s_%s' %(model.org,model.id), promseq(cdsid(model)))

# history -g -f hist # ipython

# ...map orthologous/homologous genes from reciprocal dicts into db...
# ...smart/greedy extension of putative real CDS's and resulting plausible true promoter coords by protein alignment across homologs?
# ...write promoters for sets of sets of orthologous/homologous genes...
# Lhcs
# ACT1's
# EFs
# metabolic pathways
# carbon concentrating?
# anything with a cool signal in the expression clustering
