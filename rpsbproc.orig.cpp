/*$Id$
*  =========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*  ===========================================================================
* 
*  Authors:  Shennan Lu
*
*  =======================================================================*/

#include <ncbi_pch.hpp>
#if defined(__BLAST_XML2__)
#include <objects/blastxml2/blastxml2__.hpp>
#else
#include <objects/blastxml/blastxml__.hpp>
#endif
#include <objects/seqloc/Na_strand.hpp>
#include <serial/objostrxml.hpp>
#include <serial/objistr.hpp>
#include <serial/serial.hpp>
#include <corelib/ncbireg.hpp>
#include <algorithm>
#include <vector>
#include <string>
#include <list>
#include <stack>
#include <map>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>

USING_NCBI_SCOPE;
using namespace objects;

/**********************************************************************
*	Constants defined here
***********************************************************************/

const string k_strEmptyString("");
	
const char * const DATAPATH = "datapath";
const char * const CDDIDS = "cdd";
const char * const CDTRACKINFO = "cdt";
const char * const CLUSTERLINKS = "clst";
const char * const FEATURES = "feats";
const char * const GENERIC_FEATURES = "genfeats";
const char * const SPECIFICTHRESHOLDS = "spthr";

const int COORDSBASE = 1;

const char DELIMIT = '\t';
const char COORDELIMIT = ',';

const char * const HITTYPE_SPECIFIC = "Specific";
const char * const HITTYPE_NONSPECIFIC = "Non-specific";
const char * const HITTYPE_CLUSTER = "Superfamily";
const char * const HITTYPE_MULTIDOM = "Multidom";
const char * const ANNOTTYPE_SPECIFIC = "Specific";
const char * const ANNOTTYPE_GENERIC = "Generic";
const char * const QUERY_TYPE_PEPTIDE = "Peptide";
const char * const QUERY_TYPE_NUCLEOTIDE = "Nucleotide";

const char * const PROGRAM_TITLE = "Post-RPSBLAST Processing Utility v0.11";

const char * const DATASTART = "DATA";
const char * const DATAEND = "ENDDATA";

const char * const SESSIONSTART = "SESSION";
const char * const SESSIONEND = "ENDSESSION";

const char * const QUERYSTART = "QUERY";
const char * const QUERYEND = "ENDQUERY";

const char * const DOMSTART = "DOMAINS";
const char * const DOMEND = "ENDDOMAINS";

const char * const FAMSTART = "SUPERFAMILIES";
const char * const FAMEND = "ENDSUPERFAMILIES";

const char * const FEATSTART = "SITES";
const char * const FEATEND = "ENDSITES";

const char * const MOTIFSTART = "MOTIFS";
const char * const MOTIFEND = "ENDMOTIFS";

const char * const CONFIGFILE = "rpsbproc.ini";
// -- standard file names
const char * const CDIDFILE = "cddid.tbl";
const char * const CDTRACKFILE = "cdtrack.txt";
const char * const CLSTLINKFILE = "family_superfamily_links.txt";
const char * const SPFEATFILE = "cddannot.dat";
const char * const GENFEATFILE = "cddannot_generic.dat";
const char * const MINBSCOREFILE = "bitscore_specific.txt";

const int OVERLAPLEADING = 500000;
/**********************************************************************
*	Commandline switches defined here
***********************************************************************/
struct TCMLSwitches
{
	// -- called to actually process the switch
	// -- first parameter will be the index of the switch. Program should
	// -- give an enum type to name indice of all switches for better sense.
	// -- if invalid (ie undefined) switches are found, first argument will be
	// -- m_iTotalSwitches. If a parameter is not a switch, the first argument
	// -- will be -1.
	typedef void lpfnSwitchProcessor(int, const char * );
	struct TCMLSwitch
	{
		char m_cToken;	//the switch token
		bool m_bUseParam;	//have a parameter
		
		TCMLSwitch(char tkn, bool useparam):
			m_cToken(tkn), m_bUseParam(useparam)
		{};
	};
	
	static const int NONSWITCH = -1;
	
	const TCMLSwitch *m_dimSWITCHES;
	int m_iTotalSwitches;
	const char * m_dimSwitchChar;	//char that identifies a switch, such as / or -. Must terminate with '\0'
	TCMLSwitches(const TCMLSwitch *pSwDefs, int ttl, const char * sc = "-"): m_dimSWITCHES(pSwDefs), m_iTotalSwitches(ttl), m_dimSwitchChar(sc) {};
	bool IsSwitchChar(char c) const;
	// -- return index of the switch. If not found, return m_iTotalSwitches
	// -- as an indication of invalid switch. bUseParam set to false.
	int FindSwitch (char cToken, bool &bUseParam) const;
	void IterParameters(int argc, char * argv[], lpfnSwitchProcessor proc) const;
		
	
};


bool TCMLSwitches::IsSwitchChar(char c) const
{
	for (const char *p = m_dimSwitchChar; 0 != *p; ++p)
		if (*p == c) return true;
	return false;
}

// -- return index of the switch. If not found, return m_iTotalSwitches
// -- as an indication of invalid switch. bUseParam set to false.
int TCMLSwitches::FindSwitch (char cToken, bool &bUseParam) const
{
	for (int i = 0; i < m_iTotalSwitches; ++i)
	{
		if (m_dimSWITCHES[i].m_cToken == cToken)
		{
			bUseParam = m_dimSWITCHES[i].m_bUseParam;
			return i;
		}
	}
	
	bUseParam = false;
	return m_iTotalSwitches;	//invalid switch
}

void TCMLSwitches::IterParameters(int argc, char * argv[], lpfnSwitchProcessor proc) const
{
	int cmdStatus = NONSWITCH;
	bool bReqParam = false;
	
	for (int i = 1; i < argc; ++i)
	{
		size_t argCharCnt = strlen(argv[i]);
		size_t argCharIdx = 0;
		
		if (argCharIdx < argCharCnt)	//has content
		{
			if (IsSwitchChar(argv[i][0]))	//is switch
			{
				++argCharIdx;
				if (argCharIdx < argCharCnt)
				{
					// -- new switch, end of last switch no matter what
					if (NONSWITCH != cmdStatus)
					{
						proc(cmdStatus, NULL);
						cmdStatus = NONSWITCH;
					}
					// -- find it
					cmdStatus = FindSwitch(argv[i][argCharIdx], bReqParam);

					
					if (bReqParam)
					{
						++argCharIdx;
						if (argCharIdx < argCharCnt)	//has immediately followed parameter
						{
							proc(cmdStatus, argv[i] + argCharIdx);
							cmdStatus = NONSWITCH;
						}
						else	//wait for next parameter as 
							continue;
					}
					else	//no parameter, check any other
					{
						// call
						if (cmdStatus == m_iTotalSwitches)	//invalid
							proc(cmdStatus, argv[i] + argCharIdx);
						else
							proc(cmdStatus, NULL);
						cmdStatus = NONSWITCH;
						++argCharIdx;
						while (argCharIdx < argCharCnt)
						{
							cmdStatus = FindSwitch(argv[i][argCharIdx], bReqParam);
							if (bReqParam)	//requires parameter, the rest of string to be parameter
							{
								++argCharIdx;
								if (argCharIdx < argCharCnt)
								{
									proc(cmdStatus, argv[i] + argCharIdx);
									cmdStatus = NONSWITCH;
								}
								else	//wait for next commandline
									break;
							}
							else	//no parameters needed
							{
								if (cmdStatus == m_iTotalSwitches)	//invalid
									proc(cmdStatus, argv[i] + argCharIdx);
								else
									proc(cmdStatus, NULL);
								cmdStatus = NONSWITCH;
								++argCharIdx;
							}
						}
					}
				}
				else	//no switch, just a '-' or '/', treat as non-switch
				{
					proc(cmdStatus, argv[i]);
					cmdStatus = NONSWITCH;
				}
			}
			else	//not a switch
			{
				proc(cmdStatus, argv[i]);	//call as a parameter for last switch (need a parameter)
				cmdStatus = NONSWITCH;
			}
		}
	}
}


/**********************************************************************
*	Commandline switches - numeric code
***********************************************************************/
struct TRpsbProcCmds: public TCMLSwitches
{
	enum EIndex
	{
		eCfgFile = 0,
		eInFile = eCfgFile + 1,	//source data stream (blast results in xml format), default to stdin
		eOutFile = eInFile + 1,	//destination (data output), default to stdout
		eEVCutoff = eOutFile + 1,	//filter with evalue
		eMode = eEVCutoff + 1,	//concise (default), standard and full
		eTData = eMode + 1,
		eData = eTData + 1,	//path to data files
		eFams = eData + 1,
		eQuiet = eFams + 1,	//added as user request, append all domains superfamily information.
		eHelp = eQuiet + 1,	//help
		eInvalid = eHelp + 1
	};
	
	static TCMLSwitch m_dimRpsbProcSwitches[eInvalid];
	
	TRpsbProcCmds(void): TCMLSwitches(m_dimRpsbProcSwitches, eInvalid) {};
	
};

TCMLSwitches::TCMLSwitch TRpsbProcCmds::m_dimRpsbProcSwitches[] = 
{
	TCMLSwitch('c', true),
	TCMLSwitch('i', true),
	TCMLSwitch('o', true),
	TCMLSwitch('e', true),
	TCMLSwitch('m', true),
	TCMLSwitch('t', true),
	TCMLSwitch('d', true),	//datapath, requested by user siewyit@ebi.ac.uk 
	TCMLSwitch('f', true),	//family
	TCMLSwitch('q', false),
	TCMLSwitch('h', false)
};

/**********************************************************************
*	Enum-Lit Type processor
***********************************************************************/
template<class EnumLitType>
int GetIdx(const string& strLit)
{
	for (int i = EnumLitType::eEnumStart; i < EnumLitType::eEnumEnd; ++i)
	{
		if (strLit == EnumLitType::dimLits[i - EnumLitType::eEnumStart]) return i;
	}
	
	return EnumLitType::eDefault;
}

/**********************************************************************
*	Result redundent level
***********************************************************************/
struct TDataModes
{
	enum EIndex
	{
		eEnumStart = 0,
		e_rep = eEnumStart,
		e_std = e_rep + 1,
		e_full = e_std + 1,
		eEnumEnd = e_full + 1
	};
	
	static const EIndex eDefault = e_rep;
	static const char* dimLits[eEnumEnd - eEnumStart];
	static const char* dimDisplay[eEnumEnd - eEnumStart];
};

const char * TDataModes::dimLits[] = {"rep", "std", "full"};
const char * TDataModes::dimDisplay[] = {"Concise", "Standard", "Full"};
	
/**********************************************************************
*	Target data
***********************************************************************/
struct TTargetData
{
	enum EIndex
	{
		eEnumStart = 0,
		e_doms = eEnumStart,
		e_feats = e_doms + 1,
		e_both = e_feats + 1,
		eEnumEnd = e_both + 1
	};
	
	static const EIndex eDefault = e_both;
	static const char* dimLits[eEnumEnd - eEnumStart];
	static const char* dimDisplay[eEnumEnd - eEnumStart];
};

const char * TTargetData::dimLits[] = {"doms", "feats", "both"};
const char * TTargetData::dimDisplay[] = {"Domain hits", "Site annotations", "Domain hits and site annotations"};


/**********************************************************************
*	Helper classes -- sequence Segment
***********************************************************************/
enum EResidueExtension
{
	eFullLeft = -2,
	eHalfLeft = -1,
	eAccurate = 0,
	eHalfRight = 1,
	eFullRight = 2,
	eUseDefaultExt = 255
	
};



struct TSeg
{
	int from;
	int to;
	
	EResidueExtension lext;
	EResidueExtension rext;
	
	// -- give up tracking conversion, just keep the original in the segment.
	int ori_from;
	
	bool operator == (const TSeg& other) const;
	
	TSeg(int f = 0, int t = 0, EResidueExtension le = eAccurate, EResidueExtension re = eAccurate): from(f), to(t), lext(le), rext(re), ori_from(f) {};
	bool IsValid(void) const;
	
	bool LeftTo(const TSeg& rSrcSeg) const {return (to < rSrcSeg.from - 1);} 
	bool RightTo(const TSeg& rSrcSeg) const {return (from > rSrcSeg.to + 1);}
	bool LeftTouch(const TSeg& rSrcSeg) const {return (to == rSrcSeg.from - 1);}
	bool RightTouch(const TSeg& rSrcSeg) const {return (from == rSrcSeg.to + 1);}
	bool MoreLeft(const TSeg& rSrcSeg) const {return (from < rSrcSeg.from);}
	bool MoreRight(const TSeg& rSrcSeg) const {return (to > rSrcSeg.to);}
	bool Overlap(const TSeg& rSrcSeg) const {return ((from >= rSrcSeg.from && from <= rSrcSeg.to) || (rSrcSeg.from >= from && rSrcSeg.from <= to));}
	bool Touch(const TSeg& rSrcSeg) const {return (LeftTouch(rSrcSeg) || RightTouch(rSrcSeg));}
	
	
	//void SaveToCache(CDataCache &dc) const;
	//void RestoreFromCache(CDataCache &dc);

};

// -- sort from left most to right
struct TSegSortLeft
{
	bool operator () (const TSeg* s1, const TSeg* s2);
};


// -- sort from long to short segs, so short segs will be drawn later than long segs

struct TSegSortLength
{
	bool operator () (const TSeg* s1, const TSeg* s2);
};



typedef list<TSeg> TSegs;

class CSegSet
{
	friend struct TOflAlignInfo;
public:
	struct TResiduePos
	{
		int curr, ori;
		TResiduePos(int c = 0, int o = 0): curr(c), ori(o) {};
	};
	
	CSegSet(void): m_iFactor(1), m_lstContainer(), m_ulGapThreshold(-1) {};
	CSegSet(const vector<int> &residues);
	CSegSet(const vector<int>& starts, const vector<unsigned int>& lens);	//always set ofs to zero! so must use slave coordinates.
	virtual ~CSegSet(void) {};
	
	void SetData(const vector<int>& starts, const vector<unsigned int>& lens);	//always set ofs to zero! so must use slave coordinates.
	void SetData(const vector<int> &residues);

	bool operator == (const CSegSet& other) const;
	
	// -- status
	bool IsEmpty(void) const {return m_lstContainer.empty();}
	int GetTransFactor(void) const {return m_iFactor;}
	
	// -- manipulate. 
	// -- any operation, the ori_from and ori_to are calculated based on target segment. src segs ori information are discarded.
	void AddSeg(int f, int t, EResidueExtension le = eAccurate, EResidueExtension re = eAccurate);	//  
	void AddSeg(const TSeg& seg);
	
	void Clear(void) {m_lstContainer.clear();}
	
	void Merge(const CSegSet& src);
	void Cross(const CSegSet& src);
	void Clip(const CSegSet& src);
	void Inv(int from, int to);	//inverse against a total range

	int GetLeft(void) const;
	int GetRight(void) const;
	int GetTotalResidues(void) const;

	void GetOverall(TSeg &target, int &ori_to) const;	//ori-to is calculated from the last segment factor is 1 or 3 -- if Pr2na
	int GetOriTo(TSegs::const_iterator citer, int pos = -1) const;
	int GetOriTo(void) const;
	const TSegs& GetSegs(void) const {return m_lstContainer;}
	
	// -- gapThreshold is AA residue counts. Will automatically convert to NA (if applicable) counts in the result
	void GetGaps(CSegSet &dst) const;
	
	void GetTranslatedPosMap(size_t aaSeqLen, std::vector<TResiduePos> &dst) const;
	void GetSimplePosMap(vector<TResiduePos> &dst) const;
	void GetOriSimplePosMap(vector<TResiduePos> &dst) const;
	
	
	virtual int GetCompleteSize(void) const {return -1;}
	
private:
	// -- the actual container
	int m_iFactor;	//factor: when map protein to na, it is 3, otherwise it is 1. sign denotes the direction
	TSegs m_lstContainer;
public:
	unsigned int m_ulGapThreshold;

};

bool TSeg::IsValid(void) const
{
	return (from >= 0 && to >= from);
}

bool TSeg::operator == (const TSeg& other) const
{
	return from == other.from && to == other.to;
}

bool TSegSortLeft::operator () (const TSeg* s1, const TSeg* s2)
{
	return s1->from < s2->from;
}

bool TSegSortLength::operator () (const TSeg* s1, const TSeg* s2)
{
	return s1->to - s1->from > s2->to - s2->from;
}

CSegSet::CSegSet(const vector<int> &residues)
{
	SetData(residues);
}

CSegSet::CSegSet(const vector<int> &starts, const vector<unsigned int> &lens): m_iFactor(1), m_lstContainer(), m_ulGapThreshold(-1)
{
	// -- do not clear segset first. we may need to add segs
	SetData(starts, lens);
}

bool CSegSet::operator == (const CSegSet& other) const
{
	return m_lstContainer == other.m_lstContainer;
}

void CSegSet::SetData(const vector<int>& starts, const vector<unsigned int>& lens)
{
	m_lstContainer.clear();
	for (size_t i = 0; i < lens.size(); ++i)
		AddSeg(starts[i], starts[i] + lens[i] - 1, eAccurate, eAccurate);
}

void CSegSet::SetData(const vector<int> &residues)
{
	m_lstContainer.clear();
	
	if (residues.empty()) return;
	
	vector<int> sorted(residues);
	sort(sorted.begin(), sorted.end());
	
	TSeg seg;
	
	seg.from = seg.to = seg.ori_from = sorted[0];
	
	size_t i = 0, iEnd = sorted.size();
	
	while (i < iEnd)
	{
		if (sorted[i] > seg.to + 1)
		{
			AddSeg(seg);
			seg.from = seg.to = seg.ori_from = sorted[i];
		}
		else if (sorted[i] > seg.to)
		{
			seg.to = sorted[i];
		}
		
		++i;
	}
	AddSeg(seg);
}


void CSegSet::AddSeg(int f, int t, EResidueExtension le, EResidueExtension re)
{
	AddSeg(TSeg(f, t, le, re));
}

// -- will discard seg original info
void CSegSet::AddSeg(const TSeg& seg)
{
	if (!seg.IsValid()) return;
	TSegs::iterator curr = m_lstContainer.begin(), dstEnd = m_lstContainer.end();
	
	while (dstEnd != curr && curr -> to < seg.from - 1) ++curr;
	if (curr == m_lstContainer.end())	//add to end
	{
		m_lstContainer.insert(dstEnd, seg);
		return;
	}
	else if (curr->from > seg.to + 1)	// gap big enough
	{
		m_lstContainer.insert(curr, seg);
		return;
	}
	else	//merge the seg in
	{
		if (curr->from > seg.from)
		{
			curr->ori_from -= (curr->from - seg.from) / m_iFactor;
			curr->from = seg.from;	//extended left
			curr->lext = seg.lext;
			
		}
		if (curr->to < seg.to)
		{
			curr->to = seg.to;	//extended right
			curr->rext = seg.rext;
		}
		
		TSegs::iterator nxt = curr;
		
		++nxt;
		while (nxt != dstEnd && nxt->from <= curr->to + 1)
		{
			if (curr->to < nxt->to)
			{
				curr->to = nxt->to;	//merge next
				curr->rext = nxt->rext;
			}
			m_lstContainer.erase(nxt);
			nxt = curr;
			++nxt;
		}
	}
}



void CSegSet::Merge(const CSegSet& src)
{
	TSegs::iterator curr = m_lstContainer.begin();
	TSegs::const_iterator scur = src.m_lstContainer.begin();
	
	while (curr!= m_lstContainer.end() && scur != src.m_lstContainer.end())
	{
		if (curr->to < scur->from - 1) ++curr;
		else if (scur->to < curr->from - 1)
		{
			m_lstContainer.insert(curr, *scur);
			++scur;
		}
		else
		{
			// -- merge this src seg
			if (curr->from > scur->from)
			{
				curr->ori_from -= (curr->from - scur->from) / m_iFactor;
				curr->from = scur->from;
				curr->lext = scur->lext;
			}
			if (curr->to < scur->to)
			{
				curr->to = scur->to;
				curr->rext = scur->rext;
			}
			// -- to next src seg
			++scur;
			
			// maintain this set			
			TSegs::iterator nxt = curr;
			++nxt;
			while (nxt != m_lstContainer.end() && nxt->from <= curr->to + 1)
			{
				if (curr->to < nxt->to)
				{
					curr->to = nxt->to;	//merge next
					curr->rext = nxt->rext;
				}
				m_lstContainer.erase(nxt);
				nxt = curr;
				++nxt;
			}
		}
	}
	// -- add the rest segs in
	
	while (scur != src.m_lstContainer.end())
		m_lstContainer.push_back(*scur++);
}

void CSegSet::Cross(const CSegSet& src)
{
	TSegs::iterator curr = m_lstContainer.begin();
	TSegs::const_iterator scur = src.m_lstContainer.begin();
	
	while (curr!= m_lstContainer.end() && scur != src.m_lstContainer.end())
	{
		if (curr->to < scur->from)
		{
			TSegs::iterator del = curr;
			++curr;
			m_lstContainer.erase(del);
		} 
		else if (scur->to < curr->from)
		{
			++scur;
		}
		else
		{
			if (curr->from < scur->from)
			{
				curr->ori_from += (scur->from - curr->from) / m_iFactor;
				curr->from = scur->from;
				curr->lext = scur->lext;
			}
			if (curr->to > scur->to + 1)
			{
				TSeg temp(curr->from, scur->to, curr->lext, scur->rext);
				temp.ori_from = curr->ori_from;
				m_lstContainer.insert(curr, temp);
				
				curr->ori_from += (scur->to + 2 - curr->from) / m_iFactor;
				curr->from = scur->to + 2;
				++scur;
			}
			else
			{
				if (curr->to > scur->to)
				{
					curr->to = scur->to;
					curr->rext = scur->rext;
					++scur;
				}
				++curr;
			}
		}
	}
	m_lstContainer.erase(curr, m_lstContainer.end());
}

void CSegSet::Clip(const CSegSet& src)
{
	TSegs::iterator curr = m_lstContainer.begin();
	TSegs::const_iterator scur = src.m_lstContainer.begin();
	
	while (curr!= m_lstContainer.end() && scur != src.m_lstContainer.end())
	{
		if (curr->to < scur->from)
		{
			++curr;
		} 
		else if (scur->to < curr->from)
		{
			++scur;
		}
		else
		{
			if (curr->from < scur->from)
			{
				TSeg temp(curr->from, scur->from - 1, curr->lext, eFullRight);
				temp.ori_from = curr->ori_from;
				m_lstContainer.insert(curr, temp);
			}
			if (curr->to <= scur->to)
			{
				TSegs::iterator del = curr;
				++curr;
				m_lstContainer.erase(del);
			}
			else
			{
				curr->ori_from += (scur->to + 1 - curr->from) / m_iFactor;
				curr->from = scur->to + 1;
				curr->lext = eFullLeft;
				++scur;
			}
		}
	}
}

void CSegSet::Inv(int from, int to)	//inverse against a total range, discard original coordinates
{
	TSegs::iterator hdr = m_lstContainer.begin(), curr = hdr, dstEnd = m_lstContainer.end();
	
	while (curr != dstEnd && curr->to < from)
		++curr;
	m_lstContainer.erase(hdr, curr);
	
	int start = from;
	while (curr != dstEnd && curr->from < to)
	{
		if (curr->from > start)	//insert new
		{
			TSeg temp(start, curr->from - 1);
			m_lstContainer.insert(curr, temp);
		}
		start = curr->to + 1;
		TSegs::iterator del = curr;
		++curr;
		m_lstContainer.erase(del);
	}
	
	// -- assert()
	if (start <= to)
		m_lstContainer.insert(curr, TSeg(start, to));
	// -- assert(curr->from > end);
	m_lstContainer.erase(curr, dstEnd);
}


int CSegSet::GetLeft(void) const
{
	return IsEmpty() ? -1 : m_lstContainer.begin()->from;
}

int CSegSet::GetRight(void) const
{
	return IsEmpty() ? -1 : m_lstContainer.rbegin()->to;
}

int CSegSet::GetTotalResidues(void) const
{
	int iResult = 0;
	for (TSegs::const_iterator iter = m_lstContainer.begin(); iter != m_lstContainer.end(); ++iter)
	{
		iResult += iter->to - iter->from + 1;
	}
	return iResult/abs(m_iFactor);
}


void CSegSet::GetOverall(TSeg &target, int &ori_to) const
{
	if (IsEmpty())
	{
		target.from = target.to = target.ori_from = -1;
	}
	else
	{
		TSegs::const_iterator iterFirst = m_lstContainer.begin();
		TSegs::const_reverse_iterator iterLast = m_lstContainer.rbegin();
		target.from = iterFirst->from;
		target.ori_from = iterFirst->ori_from;
		target.lext = iterFirst->lext;
		target.to = iterLast->to;
		target.rext = iterLast->rext;
		ori_to = iterLast->ori_from + (iterLast->to - iterLast->from) / m_iFactor;
	}
}

int CSegSet::GetOriTo(void) const
{
	if (m_lstContainer.empty()) return -1;
	TSegs::const_reverse_iterator riter = m_lstContainer.rbegin();
		
	return riter->ori_from + (riter->to - riter->from) / m_iFactor;
}

int CSegSet::GetOriTo(TSegs::const_iterator citer, int pos) const
{
	//return citer->ori_from + (citer->to - citer->from) / m_iFactor;
	if (pos < 0) pos = citer->to;
	return citer->ori_from + (pos - citer->from) / m_iFactor;
}

void CSegSet::GetGaps(CSegSet &dst) const
{
	dst.Clear();
	
	if (IsEmpty() || ((size_t)-1 == m_ulGapThreshold)) return;
	int fac = abs(m_iFactor);
	unsigned int gapThreshold = fac * m_ulGapThreshold;
	if (3 == fac && gapThreshold > 2) gapThreshold -= 2;
	
	TSegs::const_iterator iter = m_lstContainer.begin(), iter2 = iter, iterEnd = m_lstContainer.end();
	
	++iter2;
	
	
	while (iterEnd != iter2)
	{
		int f = iter->to + 1, t = iter2->from - 1;
		if (t - f + 1 > (int)gapThreshold)	//consider as gap
			dst.AddSeg(f, t, eFullLeft, eFullRight);
			
		++iter;
		++iter2;
	}
}

void CSegSet::GetTranslatedPosMap(size_t aaSeqLen, vector<CSegSet::TResiduePos> &dst) const
{
	if (m_iFactor > 0)	//positive, nothing to worry about
	{
		for (TSegs::const_iterator iter = m_lstContainer.begin(), iterEnd = m_lstContainer.end(); iter != iterEnd; ++iter)
		{
			for (TSignedSeqPos c = iter->from, inc = 0; c <= iter->to; c += m_iFactor, ++inc)
			{
				dst.push_back(TResiduePos(c / m_iFactor, iter->ori_from + inc));
			}
		}
	}
	else
	{
		for (TSegs::const_iterator iter = m_lstContainer.begin(), iterEnd = m_lstContainer.end(); iter != iterEnd; ++iter)
		{
			for (TSignedSeqPos c = iter->to, inc = 0; c > iter->from; c += m_iFactor, ++inc)
			{
				dst.push_back(TResiduePos(aaSeqLen + c / m_iFactor, iter->ori_from + inc));
			}
		}
	}
}

void CSegSet::GetSimplePosMap(vector<CSegSet::TResiduePos> &dst) const
{
	dst.clear();
//	int sc = m_iFactor > 0 ? m_iFactor : -m_iFactor;
	for (TSegs::const_iterator iter = m_lstContainer.begin(), iterEnd = m_lstContainer.end(); iter != iterEnd; ++iter)
	{
		for (int c = iter->from; c <= iter->to; ++c)
		{
			dst.push_back(TResiduePos(c, iter->ori_from + (c - iter->from) / m_iFactor));
		}
	}
}


void CSegSet::GetOriSimplePosMap(vector<CSegSet::TResiduePos> &dst) const
{
	dst.clear();
	int fact = abs(m_iFactor);
	
	for (TSegs::const_iterator iter = m_lstContainer.begin(), iterEnd = m_lstContainer.end(); iter != iterEnd; ++iter)
	{
		for (int c = iter->from; c <= iter->to; c += fact)
		{
			dst.push_back(TResiduePos(c, iter->ori_from + (c - iter->from) / m_iFactor));
		}
	}
}


/**********************************************************************
*	Helper class -- CProSite
***********************************************************************/
class CProSite
{
public:
	enum EParseError
	{
		eNoError,
		eInvalidChar,
		eSyntaxError,
		eUnexpectedEnd
	};

	struct TMatched
	{
		size_t start, end;
		TMatched(size_t s = std::string::npos, size_t e = std::string::npos): start(s), end(e) {};
	};
	CProSite(void):m_vecCompiledPattern() {};
	EParseError Parse(const std::string &expr, size_t &errorPos);
	size_t Match(const std::string &text, size_t tlen, size_t start_pos) const;
	void Search(const std::string &text, std::vector<TMatched> &result) const;
	
	// -- return a vector of minimal match length, with 'X' as generic, 'A' as alternative and S as strict.
	void GetMinimalXMap(std::string &minMap) const;
	
private:
	struct TPatternPos
	{
		static const unsigned int NEGATIVE_FILTER = 0x1;
		static const unsigned int LAZY_MATCH = 0x2;
		unsigned int m_uiFlags;	//positive or negative select
		std::string m_strAltChars;
		size_t m_ulMinCount, m_ulMaxCount;
		size_t (TPatternPos::*m_pfnMatch)(const std::string &text, size_t tlen, size_t pos) const;
		
		TPatternPos(void);
		//return <0 as match fail. return 0 as "just match" without flex. return >0 as how many flex positions
		int Match(const std::string &text, size_t tlen, std::vector<size_t> &rMatchRec, size_t start_pos = std::string::npos) const;
		size_t x_NegMatch(const std::string &text, size_t tlen, size_t pos) const;
		size_t x_PosMatch(const std::string &text, size_t tlen, size_t pos) const;
			
		void x_NormalizeAltChars(void);
	};
	


	struct x_TMatchRec
	{
		std::vector<size_t> flex;
		std::vector<TPatternPos> :: const_iterator iterPos;
		
		x_TMatchRec(std::vector<TPatternPos> :: const_iterator p): flex(), iterPos(p) {};
	};
	
	std::vector<TPatternPos> m_vecCompiledPattern;
	
};


CProSite::TPatternPos::TPatternPos(void):
	m_uiFlags(0), m_strAltChars(k_strEmptyString), m_ulMinCount(1), m_ulMaxCount(1), m_pfnMatch(&TPatternPos::x_PosMatch)
{};



void CProSite::Search(const std::string &text, vector<CProSite::TMatched> &result) const
{
	result.clear();
	size_t tlen = text.size();
	for (size_t i = 0; i < tlen; ++i)
	{
		size_t mpos = Match(text, tlen, i);
		if (string::npos != mpos)	//success
			result.push_back(TMatched(i, mpos));
	}
}

 
CProSite::EParseError CProSite::Parse(const string &expr, size_t &errorPos)
{
	if (expr.empty()) return eUnexpectedEnd;
	
	enum __EParserStatus
	{
		eParserReady = 0,
		ePlainText,	//collecting plain text
		eAlternates,	//in [] collecting alternative
		eAltRedund,
		eNegations,	//in {} collecting alternative
		eNegRedund,
		eRepeats,	//in () collecting repeats
		eRange	//only valid for X, in (a,b) collecting 
	} status = eParserReady;
	
	EParseError eErrorType = eNoError;	

	
	string::const_iterator iterChar = expr.begin(), iterCharEnd = expr.end();
	
	m_vecCompiledPattern.clear();
	TPatternPos __dummy;
	vector<TPatternPos> :: iterator iterCurr = m_vecCompiledPattern.end();
	
	string lit(k_strEmptyString);
	bool bEscaped = false;	//if '\' is detected
	bool bInfLit = false;
	
	while (iterCharEnd != iterChar)
	{
		char c = *iterChar;
		if ('a' <= c && 'z' >= c) c -= 0x20;	//turn to capital letters
		switch (c)
		{
		case '*':	//stop codon
			if (!bEscaped)	//not escaped, as equivalent to (0, inf)
			{
				switch (status)
				{
				case ePlainText:
					iterCurr->m_ulMinCount = 0;
					iterCurr->m_ulMaxCount = -1;	//the biggest possible
					if (iterCurr->m_strAltChars == "X") iterCurr->m_uiFlags |= TPatternPos::LAZY_MATCH;
					break;
				case eAlternates:	//always treat as a literal * (stop codon)
				case eNegations:
					iterCurr->m_strAltChars.push_back(c);
					break;
				case eRange:	//accepting max
				case eRepeats:
					if (lit.empty() && !bInfLit) 
					{
						bInfLit = true;
						break;
					}
				case eAltRedund:
				case eNegRedund:
					break;
				default:
					eErrorType = eSyntaxError;
					goto labelErrorRet;
				}
				break;
			}
			
		case 'A':
		case 'C':
		case 'D':
		case 'E':
		case 'F':
		case 'G':
		case 'H':
		case 'I':
		case 'K':
		case 'L':
		case 'M':
		case 'N':
		case 'P':
		case 'Q':
		case 'R':
		case 'S':
		case 'T':
		case 'U':	//newly added selenocysteine
		case 'V':
		case 'W':
		case 'Y':
		
			switch (status)
			{
			case eParserReady:	//No current
				status = ePlainText;
			case ePlainText:	//current is valid
				iterCurr = m_vecCompiledPattern.insert(m_vecCompiledPattern.end(), __dummy);
			case eAlternates:	//single char alternative, assume: iterCurrAltPattern is undefined
			case eNegations:
				iterCurr->m_strAltChars.push_back(c);
				break;
			default:
				eErrorType = eSyntaxError;
				goto labelErrorRet;
			}
			bEscaped = false;
			break;
		case 'X':	//special: wildcard, no alternates or negates allowed.
			switch (status)
			{
			case eParserReady:	//No current
				status = ePlainText;
			case ePlainText:	//current is valid
				iterCurr = m_vecCompiledPattern.insert(m_vecCompiledPattern.end(), __dummy);
				iterCurr->m_strAltChars.push_back(c);
				//iterCurr->m_uiFlags |= TPatternPos::LAZY_MATCH;
				break;
			case eAlternates:	//if alternates, and 'X' mean everything canbe accepted.
			case eNegations:
				iterCurr->m_strAltChars.push_back(c);
				break;
			default:
				eErrorType = eSyntaxError;
				goto labelErrorRet;
			}
			bEscaped = false;
			break;
		case '0':
		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
		case '6':
		case '7':
		case '8':
		case '9':
			switch (status)
			{
			case eRepeats:
			case eRange:
				if (!bInfLit)
				{
					lit.push_back(c);
					break;
				}
			
			default:
				eErrorType = eSyntaxError;
				goto labelErrorRet;
				
			case eAltRedund:
			case eNegRedund:
				break;	//ignore it
			}
			bEscaped = false;
			break;
		case '(':
			switch (status)
			{
			case ePlainText:
				status = eRepeats;	//change to repeats. assert(cTargetChar is defined)
				lit.clear();
				break;
			// -- added 02/12/2016 to tolerate redunence in alternations
			case eAlternates:
				status = eAltRedund;
				break;
			case eNegations:
				status = eNegRedund;
				break;
			// -- End added 02/12/2016 
			default:
				eErrorType = eSyntaxError;
				goto labelErrorRet;
			}
			break;
		case ')':
			switch (status)
			{
			case eRepeats:	//end of repeat
				if (!lit.empty())	//has input
				{
					iterCurr->m_ulMinCount = iterCurr->m_ulMaxCount = (size_t)atol(lit.c_str());
				}
				else if (bInfLit)
				{
					iterCurr->m_ulMinCount = 0;
					iterCurr->m_ulMaxCount = -1;	//the biggest possible
					if (iterCurr->m_strAltChars == "X") iterCurr->m_uiFlags |= TPatternPos::LAZY_MATCH;	//infinit X will have lazy match
					bInfLit = false;
				}
				else
				{
					eErrorType = eSyntaxError;
					goto labelErrorRet;
				}
				status = ePlainText;

				break;
			case eRange:	//end of range, assume iterGroup already renewed by the ',' token
				if (!lit.empty())	//has input
				{
					iterCurr->m_ulMaxCount = (size_t)atol(lit.c_str());
					// -- check for reversed
					if (iterCurr->m_ulMaxCount < iterCurr->m_ulMinCount)
					{
						iterCurr->m_ulMaxCount ^= iterCurr->m_ulMinCount;
						iterCurr->m_ulMinCount ^= iterCurr->m_ulMaxCount;
						iterCurr->m_ulMaxCount ^= iterCurr->m_ulMinCount;
					}
				}
				else if (bInfLit)
				{
					iterCurr->m_ulMaxCount = -1;
					if (iterCurr->m_strAltChars == "X") iterCurr->m_uiFlags |= TPatternPos::LAZY_MATCH;
					bInfLit = false;
				}
				else
				{
					eErrorType = eSyntaxError;
					goto labelErrorRet;
				}
				status = ePlainText;
				break;
			case eAltRedund:
				status = eAlternates;
				break;
			case eNegRedund:
				status = eNegations;
				break;	//ignore it 
			default:
				eErrorType = eSyntaxError;
				goto labelErrorRet;
			};
			break;

		case ',':	//turn repeat into range
			if (bEscaped)
			{
				eErrorType = eInvalidChar;
				goto labelErrorRet;
			}
			switch (status)
			{
			case eRepeats:	//turn repeats into a range spec.
			
				if (!lit.empty())
				{
					iterCurr->m_ulMinCount = (size_t)atol(lit.c_str());
					lit.clear();
					status = eRange;
					break;
				}
			default:
				eErrorType = eSyntaxError;
				goto labelErrorRet;
			case eAltRedund:
			case eNegRedund:
				break;
			};
			break;
		case '\\':
			if (bEscaped)
			{
				eErrorType = eInvalidChar;
				goto labelErrorRet;
			}
			bEscaped = true;
			break;
			
		case '[':	//start variable set
			if (bEscaped)
			{
				eErrorType = eInvalidChar;
				goto labelErrorRet;
			}
			switch (status)
			{
			case ePlainText:
			case eParserReady:
				iterCurr = m_vecCompiledPattern.insert(m_vecCompiledPattern.end(), __dummy);
				//iterCurr->m_uiFlags = true;
				status = eAlternates;
				break;
			default:
				eErrorType = eSyntaxError;
				goto labelErrorRet;
			};
			break;
		case ']':	//end variable set
			if (bEscaped)
			{
				eErrorType = eInvalidChar;
				goto labelErrorRet;
			}
			switch (status)
			{
			case eAlternates:
				if (!iterCurr->m_strAltChars.empty())
				{
					iterCurr->x_NormalizeAltChars();
					status = ePlainText;
					break;
				}
				
			default:
				eErrorType = eSyntaxError;
				goto labelErrorRet;
			}
			break;
		case '{':	//negative set
			if (bEscaped)
			{
				eErrorType = eInvalidChar;
				goto labelErrorRet;
			}
			switch (status)
			{
			case ePlainText:	//start of alternative
			case eParserReady:
				iterCurr = m_vecCompiledPattern.insert(m_vecCompiledPattern.end(), __dummy);
				iterCurr->m_uiFlags |= TPatternPos::NEGATIVE_FILTER;
				iterCurr->m_pfnMatch = &TPatternPos::x_NegMatch;
				status = eNegations;
				break;
			default:
				eErrorType = eSyntaxError;
				goto labelErrorRet;
			};
			break;
		case '}':	//end negative set
			if (bEscaped)
			{
				eErrorType = eInvalidChar;
				goto labelErrorRet;
			}
			switch (status)
			{
			case eNegations:
				if (!iterCurr->m_strAltChars.empty())
				{
					iterCurr->x_NormalizeAltChars();
					status = ePlainText;
					break;
				}
			default:
				eErrorType = eSyntaxError;
				goto labelErrorRet;
			}
			break;
		case '<':	//start anchor
			if (bEscaped)
			{
				eErrorType = eInvalidChar;
				goto labelErrorRet;
			}
			switch (status)
			{
			case eAlternates:
				if (m_vecCompiledPattern.begin() != iterCurr)	//not the first 
				{
					eErrorType = eSyntaxError;
					goto labelErrorRet;
				}
			case eNegations:
				iterCurr->m_strAltChars.push_back(c);
				break;
			case eParserReady:
				iterCurr = m_vecCompiledPattern.insert(m_vecCompiledPattern.end(), __dummy);
				iterCurr->m_strAltChars.push_back(c);
				status = ePlainText;
				break;
			default:
				eErrorType = eSyntaxError;
				goto labelErrorRet;
			}
			break;
		case '>':	//end anchor
			if (bEscaped)
			{
				eErrorType = eInvalidChar;
				goto labelErrorRet;
			}
			switch (status)
			{
			case ePlainText:
			case eParserReady:
				iterCurr = m_vecCompiledPattern.insert(m_vecCompiledPattern.end(), __dummy);
				iterCurr->m_strAltChars.push_back(c);
				status = ePlainText;
				break;
			case eNegations:
				iterCurr->m_strAltChars.push_back(c);
				break;
			default:
				eErrorType = eSyntaxError;
				goto labelErrorRet;
			}
			break;
		case '-':	//ignore
		case ' ':
		case '\t':
			break;
		default:
			eErrorType = eInvalidChar;
			goto labelErrorRet;
		}
		++iterChar;
	}
	
labelErrorRet:
	errorPos = (iterChar - expr.begin());
	return eErrorType;	
}
 
size_t CProSite::Match(const string &text, size_t tlen, size_t start_pos) const
{
	if (start_pos >= tlen) return string::npos;
	stack<x_TMatchRec> stkFlexStack;
	
	x_TMatchRec curr(m_vecCompiledPattern.begin());
	
	vector<TPatternPos> :: const_iterator iterPosEnd = m_vecCompiledPattern.end();
	
	
	while (iterPosEnd != curr.iterPos)
	{
		int matched = curr.iterPos->Match(text, tlen, curr.flex, start_pos);
		if (matched < 0)	//match failed
		{
			if (stkFlexStack.empty()) return string::npos;
			curr = stkFlexStack.top();
			stkFlexStack.pop();
		}
		else	//match success
		{
			if (matched > 0) stkFlexStack.push(curr);
			start_pos = curr.flex.back();
			++curr.iterPos;
			curr.flex.clear();
		}
	}
	
	
	// -- successfully matched all PatternPosition
	return start_pos;	//should always be valid.
	
}
void CProSite::TPatternPos::x_NormalizeAltChars(void)
{
	if (!m_strAltChars.empty())
	{
		size_t idx = 0, ttl = m_strAltChars.size();
		vector<size_t> veciidx;
		veciidx.reserve(ttl);
		while (idx < ttl)
		{
			char currc = m_strAltChars[idx];
			if ('X' == currc)	//no matter others, just being the one
			{
				m_strAltChars.clear();
				m_strAltChars.push_back('X');
				break;
			}
			veciidx.clear();
			
			size_t mv = 0;
			size_t idx2 = idx + 1;
			while (idx2 < ttl)
			{
				if (m_strAltChars[idx2] == currc)
					++mv;
				else if (mv > 0)
					m_strAltChars[idx2 - mv] = m_strAltChars[idx2];
				
				++idx2;
			}
			
			ttl -= mv;
			
			while (mv > 0)
			{
				m_strAltChars.pop_back();
				--mv;
			}
			
			++idx;
		}
	}

}
// -- negative position match
size_t CProSite::TPatternPos::x_NegMatch(const string &text, size_t tlen, size_t pos) const
{
	string::const_iterator iterChar = m_strAltChars.begin(), iterCharEnd = m_strAltChars.end();
	size_t adv = 0;
	
	while (iterCharEnd != iterChar)
	{
		switch (*iterChar)
		{
		case '<':	//start anchor
			if (0 == pos)	//matched, fail it
				return string::npos;
			break;
		case '>':	//end anchor
			if (pos == tlen)	//matched
				return string::npos;
			break;
		case 'X':
			if (pos < tlen)
				return string::npos;
			break;
		default:
			if (pos < tlen)
			{
				char c = text[pos];
				adv = 1;
				if ('a' <= c && 'z' >= c) c -= 0x20;
				if (*iterChar == c)	//matched
					return string::npos;
			}
			else
				return string::npos;
		}
		
		++iterChar;
	}
	return adv;	//nothing matched, successful
}
 
// -- positive position match
size_t CProSite::TPatternPos::x_PosMatch(const string &text, size_t tlen, size_t pos) const
{
	string::const_iterator iterChar = m_strAltChars.begin(), iterCharEnd = m_strAltChars.end();
	while (iterCharEnd != iterChar)
	{
		switch (*iterChar)
		{
		case '<':	//start anchor
			if (0 == pos)	//matched
				return 0;
			break;
		case '>':	//end anchor
			if (pos == tlen)	//matched
				return 0;
			break;
		case 'X':	//match any
			if (pos < tlen)	//last, 
				return 1;
			break;
		default:
			if (pos < tlen)
			{
				char c = text[pos];
				if ('a' <= c && 'z' >= c) c -= 0x20;
				if (*iterChar == c)	//matched
					return 1;
			}
		}
		
		++iterChar;
	}
	return string::npos;
}


int CProSite::TPatternPos::Match(const string &text, size_t tlen, vector<size_t> &rMatchRec, size_t start_pos) const
{

	size_t matched = rMatchRec.size();	//assume: at least one [0] as the start position

	int retVal = 0;	//return 0: match success but no more flex, -1: match failed. >1: # of flex left
	if (m_uiFlags & LAZY_MATCH)	//do not clear rMatchRec
	{
		if (0 == matched)	//the first time
		{
			
			rMatchRec.push_back(start_pos);
			
			// -- first, must match to minimal count
			
			while (rMatchRec.size() <= m_ulMinCount)
			{
				size_t adv = (this->*m_pfnMatch)(text, tlen, start_pos);
				
				if (string::npos == adv)	//fail
				{
					return -1;
				}
				start_pos += adv;
				rMatchRec.push_back(start_pos);
			}
			if (m_ulMaxCount > m_ulMinCount) retVal = 1;
			else retVal = 0;
		}
		else if (matched <= m_ulMaxCount)	//still flex	//revisit
		{
			start_pos = rMatchRec[matched - 1];
			size_t adv = (this->*m_pfnMatch)(text, tlen, start_pos);
			if (string::npos == adv)	//successful
				retVal = -1;
			else
			{
				start_pos += adv;
				rMatchRec.push_back(start_pos);
				if (m_ulMaxCount > matched) retVal = 1;
				else retVal = 0;
			}
		}
	}
	else	//aggressive match
	{
		if (0 == matched)	//the first time
		{
			rMatchRec.push_back(start_pos);
			while ((matched = rMatchRec.size()) <= m_ulMaxCount)
			{
				size_t adv = (this->*m_pfnMatch)(text, tlen, start_pos);
				if (string::npos == adv)	//fail
					break;
				start_pos += adv;
				rMatchRec.push_back(start_pos);
			}
		}
		else 
		{
			rMatchRec.pop_back();
			matched = rMatchRec.size();
		}
		if (matched < m_ulMinCount + 1) retVal = -1;	// minimal match not reached
		else if (matched > m_ulMinCount + 1) retVal = 1;	//with flex
		else retVal = 0;
	}
	return retVal;
}

void CProSite::GetMinimalXMap(string &minMap) const
{
	minMap.clear();
	size_t ttl = m_vecCompiledPattern.size();

	if (ttl > 0)
	{
		minMap.reserve(ttl + ttl);
		for (size_t i = 0; i < ttl; ++i)
		{
			if (m_vecCompiledPattern[i].m_ulMinCount > 0)
			{
				switch (m_vecCompiledPattern[i].m_strAltChars.size())
				{
				case 0:
					break;
				case 1:
					if ('X' == m_vecCompiledPattern[i].m_strAltChars[0]) minMap.append(m_vecCompiledPattern[i].m_ulMinCount, 'X');
					else if (m_vecCompiledPattern[i].m_uiFlags & TPatternPos::NEGATIVE_FILTER) minMap.append(m_vecCompiledPattern[i].m_ulMinCount, 'A');
					else minMap.append(m_vecCompiledPattern[i].m_ulMinCount, 'S');
					break;
				default:	//>1
					minMap.append(m_vecCompiledPattern[i].m_ulMinCount, 'A');
				}
			}
		}
	}
	
}


/**********************************************************************
*	Global data
***********************************************************************/
string g_strCfgFile(k_strEmptyString);
string g_strSrcFile(k_strEmptyString);
string g_strDstFile(k_strEmptyString);
string g_strDataPath(k_strEmptyString);

double g_dEValue = 0.01;
int g_iDataMode = TDataModes::eDefault;
int g_iTargetData = TTargetData::eDefault;
bool g_bSuperfams = false;
bool g_bRunProg = true;	//by default, program will run, but if user specified -h to get help, program won't run, just display usage information
bool g_bSilent = false;	//display program info

const char * g_pProgName = NULL;



/**********************************************************************
*	NCBI Registry handler -- config file
***********************************************************************/
bool ReadConfig(const char * cfgfile, CNcbiRegistry &reg)
{
	if (!cfgfile) cfgfile = CONFIGFILE;
	bool status = false;
	
	ifstream regfs(cfgfile, ios::in | ios::binary);

	if (regfs.good())
	{
		reg.Read(regfs);
		status = true;
	}
		
	regfs.close();
	return status;
}


/**********************************************************************
*	Bio-objects in/out templates
***********************************************************************/
template<class BioObj> 
void ObjStreamIn(CNcbiIstream& is, BioObj& rBioObj, ESerialDataFormat eFormat = eSerial_AsnText)
{
	CObjectIStream *pIStream = CObjectIStream::Open(eFormat, is, eNoOwnership);
	*pIStream >> rBioObj;
	//pIStream->Read(ObjectInfo(rBioObj));
	delete pIStream;
}

template<class BioObj> 
BioObj* ObjStreamIn(CNcbiIstream& is, BioObj* pBioObj = nullptr, ESerialDataFormat eFormat = eSerial_AsnText)
{
	BioObj* pEffective = pBioObj ? pBioObj : new BioObj();
	CObjectIStream *pIStream = CObjectIStream::Open(eFormat, is, eNoOwnership);

	*pIStream >> (*pEffective);
	delete pIStream;
	return pEffective;
}

template<class BioObj> 
void ObjStreamOut(CNcbiOstream& os, const BioObj& rBioObj, ESerialDataFormat eFormat = eSerial_AsnText)
{

  if (eSerial_Xml == eFormat)	//Get rid of dtd
  {
  	CObjectOStreamXml os_xml(os, false);
  	os_xml.SetReferenceDTD(false);
  	os_xml << rBioObj;
  }
	else
	{
		CObjectOStream *pOStream = CObjectOStream::Open(eFormat, os, eNoOwnership);
		*pOStream << rBioObj;
		delete pOStream;
	}
}

template<class BioObj> 
void ObjStreamOut(CNcbiOstream& os, const BioObj* pBioObj, ESerialDataFormat eFormat = eSerial_AsnText)
{
	if (nullptr == pBioObj) return;
	
	if (eSerial_Xml == eFormat)	//Get rid of dtd
    {
    	CObjectOStreamXml os_xml(os, false);
    	os_xml.SetReferenceDTD(false);
    	os_xml << (*pBioObj);
    }
	else
	{
		CObjectOStream *pOStream = CObjectOStream::Open(eFormat, os, eNoOwnership);
		*pOStream << (*pBioObj);
		delete pOStream;
	}
}


struct TDomSrcCount
{
	enum ESrcIdx
	{
		eCDD = 0,
		ePFam = eCDD + 1,
		eTIGRFam = ePFam + 1,
		eCOG = eTIGRFam + 1,
		eSMART = eCOG + 1,
		ePRK = eSMART + 1,
		TOTALSRCS = ePRK + 1
	};
	
	//static const size_t TOTALSRCS = ePRK + 1;
	static const size_t TOTALSIGS = 11;
	static const char * DOMSRCSIGS[TOTALSIGS];
	static const int MAXCOUNTS[TOTALSRCS];
	static ESrcIdx DomAccType(const string &acxn);
	
	
	map<ESrcIdx, int> m_SrcCounter;
	// -- return true: counted. false: src already full, so not counted
	bool CountSrc(const std::string &acxn);
};

const int TDomSrcCount::MAXCOUNTS[] = {1, 1, 1, 1, 1, 1};
const char * TDomSrcCount::DOMSRCSIGS[] = {"CD", "PFAM", "TIGR", "COG", "SMART", "PRK", "CHL", "MTH", "PHA", "PLN", "PTZ"};

TDomSrcCount::ESrcIdx TDomSrcCount::DomAccType(const string &acxn)
{
	size_t iSig = 0, iChar = 0, iAccChar = 0, acclen = acxn.size();
	

	while (iSig < TOTALSIGS && 0 != DOMSRCSIGS[iSig][iChar] && iAccChar < acclen)
	{
		char accChar = acxn[iAccChar];
		if ('a' <= accChar && 'z' >= accChar) accChar -= 0x20;	//to uppercase
		if (accChar == DOMSRCSIGS[iSig][iChar])
		{
			++iAccChar;
			++iChar;
		}
		else
		{

			++iSig;
			iChar = 0;
			iAccChar = 0;
		}
	}
/**************************************************************************************
Besides Pfam, CDD, TIGRFAM, COGs, SMART, PRK, we still have these type of accessions:
"CHL", "MTH", "PHA", "PLN", "PTZ"
All of these should be grouped together with the “PRK” as “PRK” models
***************************************************************************************/

	if (iSig < ePRK)
		return static_cast<ESrcIdx> (iSig);
	else if (iSig < TOTALSIGS)
		return ePRK;
	else
		return TOTALSRCS;
}

bool TDomSrcCount::CountSrc(const std::string &acxn)
{
	ESrcIdx srctype = DomAccType(acxn);
	if (TOTALSRCS == srctype) return true;	//unknown src always enter.
	
	map<ESrcIdx, int> :: iterator iter = m_SrcCounter.emplace(srctype, 0).first;


	if (iter->second >= MAXCOUNTS[srctype]) return false;
	++iter->second;
	return true;
}




/**********************************************************************
*	Data structure and processing classes
***********************************************************************/
const int RF_SIZE = 3;	//reading frame size
const int TOTAL_RFS = 6;	//reading frame size
const char* RF_TITLES[TOTAL_RFS] = {"RF +1", "RF +2", "RF +3", "RF -1", "RF -2", "RF -3"};

inline int RF2Idx(int rf)
{
	return (rf > 0 ? rf - 1 : -rf + 2);
}

inline int Idx2RF(int idx)
{
	return (idx > 2 ? 2 - idx : idx + 1);
}

/**********************************************************************
*	Biodata structure -- domain feature (site annotation)
***********************************************************************/
struct TOflCDFeat: public CSegSet
{
	enum EFeatType	//this is arbitrary. 
	{
		eType_Other = 0,
		eType_Active = 1,
		eType_PolyPBinding = 2,
		eType_NtBinding = 3,
		eType_IonBinding = 4,
		eType_ChemBinding = 5,
		eType_PostTransMod = 6,
		eType_StructMotif = 7
	};
	
	static const int TOTAL_TYPES = eType_StructMotif + 1;
	static const char * GENERIC_SITE_TITLE;
	
	static const char * FEATTYPES[TOTAL_TYPES];
	
	static const unsigned int STRUCTURE_BASED_EVIDENCE = 0x1 << 0;
	static const unsigned int REFERENCE_BASED_EVIDENCE = 0x1 << 1;
	static const unsigned int ADDITIONAL_COMMENTS = 0x1 << 2;
	
	string m_strTitle;	//short title
	string m_strDescr;	//descrption, future use
	string m_strMotif;
	string m_strModifiedMotif;	//depends on m_iMotifuse, this is a working version for CProSite
	int m_iMotifuse;
	int m_iIndex;	//index within pssm or other sequence
	int m_iType;
	int m_iCompleteSize;
	
	unsigned int flags;
	
	
	
	virtual int GetCompleteSize(void) const {return m_iCompleteSize;}
	void SetMotifStr(const string &rMotifStr);
	
	// -- this is to checked the mapped site with motif. If rSeqData is provided, it ckecks
	// -- both the residue position and type. If rSeqData is not provided, just check for
	// -- all non-x residues are mapped. return values:
	// -- 0: success or No motif to check.
	// -- 1: Essential positions not complete
	// -- 2: Essential positions complete but residue type mismatch.
	int MotifCheck(const vector<CSegSet::TResiduePos> &rMappedRes, const string &rSeqData = k_strEmptyString) const;	//TResiduePos contains mapped and original positions, both on protein sequence
	TOflCDFeat(void);
	
};

const char * TOflCDFeat::GENERIC_SITE_TITLE = "active site";
const char * TOflCDFeat::FEATTYPES[] = {"other", "active site", "polypeptide binding site", "nucleotide binding site", "ion binding site", "chemical binding site", "posttranslational modification", "structural motif"};

TOflCDFeat::TOflCDFeat(void):
	CSegSet(), m_strTitle(k_strEmptyString), m_strDescr(k_strEmptyString), m_strMotif(k_strEmptyString), m_strModifiedMotif(k_strEmptyString), 
	m_iMotifuse(0), m_iIndex(0), m_iType(eType_Active), m_iCompleteSize(0), flags(0)
{}

void TOflCDFeat::SetMotifStr(const string &rMotifStr)
{
	m_strMotif = rMotifStr;
	
	//// -- modify. at this time, just turn all 'X' to X* so the non-essential residues are optional
	//m_strModifiedMotif.clear();
	//char cLast = 0;
	//int count = 0;
	//for (string::const_iterator iter = m_strMotif.begin(), iterEnd = m_strMotif.end(); iter != iterEnd; ++iter)
	//{
	//	char cCurr = *iter;
	//	if (cCurr >= 'a' && cCurr <= 'z') cCurr -= 0x20;
	//	
	//	if (cCurr >= 'A' && cCurr <= 'Z')	//is residue
	//	{
	//		if (cCurr == cLast)
	//		{
	//			if ('X' != cLast) ++count;
	//		}
	//		else	//cCurr != cLast
	//		{
	//			if ('X' != cLast && count > 1)
	//			{
	//				char dimBuf[16];
	//				sprintf(dimBuf, "(%d)", count);
	//				m_strModifiedMotif.append(dimBuf);
	//				
	//			}
	//			
	//			count = 0;
	//			cLast = cCurr;
	//			
	//			m_strModifiedMotif.push_back(cCurr);
	//			
	//			if ('X' == cCurr)
	//				m_strModifiedMotif.push_back('*');
	//			else
	//				count = 1;
	//				
	//		}
	//	}
	//	else	//not residue
	//	{
	//		if (count > 1 && 'X' != cLast)	//
	//		{
	//			char dimBuf[16];
	//			sprintf(dimBuf, "(%d)", count);
	//			m_strModifiedMotif.append(dimBuf);
	//		}
	//		count = 0;
	//		cLast = 0;
	//		
	//		m_strModifiedMotif.push_back(cCurr);
	//	}
	//}
}

int TOflCDFeat::MotifCheck(const vector<CSegSet::TResiduePos> &rMappedRes, const string &rSeqData) const
{
	if (!m_strMotif.empty())
	{
		size_t seqLen = rSeqData.size();
		vector<CSegSet::TResiduePos> vecOriPoses;
		GetTranslatedPosMap(seqLen, vecOriPoses);


		CProSite ps;
		size_t errPos;
		CProSite::EParseError err = ps.Parse(m_strMotif, errPos);
		if (CProSite::eNoError != err)
		{
			cerr << "Motif string parse error -- Motif = " << m_strMotif << ", error position: " << errPos;
		}
		string minMap(k_strEmptyString);
		ps.GetMinimalXMap(minMap);
		
		size_t mtfLen = minMap.size();
		size_t mappedLen = rMappedRes.size();
//		assert(mtfLen == vecOriPoses.size());
		
		for (size_t i = 0; i < mtfLen; ++i)
		{

			if (minMap[i] != 'X')	//x always match
			{
				for (size_t j = 0; j < mappedLen; ++j)
				{
					if (rMappedRes[j].ori == vecOriPoses[i].curr)	//found, means mapped
						goto labelResidueMapped;
				}
				return 1;
			labelResidueMapped:;
					
			}
		}
		
		if (!rSeqData.empty())
		{
			
			minMap.clear();	//borrow this for other use
			
			for (size_t i = 0; i < mappedLen; ++i)
				if ((size_t)rMappedRes[i].curr < seqLen)
				{

					minMap.push_back(rSeqData[rMappedRes[i].curr]);
				}

			size_t endPos = ps.Match(minMap, seqLen, 0);
			if (string::npos == endPos)
			{
				return 2;
			}
		}
	}

	return 0;


//	if (!m_strMotif.empty())
//	{
//		vector<CSegSet::TResiduePos> vecOriPoses;
//		GetSimplePosMap(vecOriPoses);
//		
//
//		CProSite ps;
//		size_t errPos;
//		CProSite::EParseError err = ps.Parse(m_strMotif, errPos);
//		if (CProSite::eNoError != err)
//		{
//			cerr << "Motif string parse error -- Motif = " << m_strMotif << ", error position: " << errPos;
//			//throw CException(CDiagCompileInfo(__FILE__, __LINE__, "int TOflCDFeat::MotifCheck(const CSegSet &mapped, const string &rSeqData) const"), NULL, CException::eUnknown, ss.str());
//		}
//		string minMap(k_strEmptyString);
//		ps.GetMinimalXMap(minMap);
//		
//		size_t mtfLen = minMap.size();
//		size_t mappedLen = rMappedRes.size();
////		assert(mtfLen == vecOriPoses.size());
//		
//		for (size_t i = 0; i < mtfLen; ++i)
//		{
//			if (minMap[i] != 'X')
//			{
//				for (size_t j = 0; j < mappedLen; ++j)
//				{
//					if (rMappedRes[j].ori == vecOriPoses[i].curr)	//found, means mapped
//						goto labelResidueMapped;
//				}
//				// --  not mapped, return failure
//				return 1;
//			labelResidueMapped:;
//					
//			}
//		}
//		
//		if (!rSeqData.empty())
//		{
//			err = ps.Parse(m_strModifiedMotif, errPos);
//			if (CProSite::eNoError != err)
//			{
//				cerr << "Motif string parse error -- Working Motif = " << m_strMotif << ", error position: " << errPos;
//			}
//			
//			minMap.clear();	//borrow this for other use
//			size_t seqLen = rSeqData.size();
//			for (size_t i = 0; i < mappedLen; ++i)
//				if ((size_t)rMappedRes[i].curr < seqLen)
//				{
//
//					minMap.push_back(rSeqData[rMappedRes[i].curr]);
//				}
//			size_t endPos = ps.Match(minMap, seqLen, 0);
//			if (string::npos == endPos)
//			{
//				return 2;
//			}
//		}
//		
//	}
//	return 0;
}

/**********************************************************************
*	Biodata structure -- Cluster (super family)
***********************************************************************/
struct TOflClusterInfo
{
	//unsigned int m_uiPSSMID;
	string m_strAccession;
	string m_strShortName;
	string m_strDefline;
	unsigned int m_uiLength;
	TOflClusterInfo(void);
	void Reset(void);
};

TOflClusterInfo::TOflClusterInfo(void):
	m_strAccession(k_strEmptyString), m_strShortName(k_strEmptyString), m_strDefline(k_strEmptyString), m_uiLength(0)
{}

void TOflClusterInfo::Reset(void)
{
	m_strAccession.clear();
	m_strShortName.clear();
	m_strDefline.clear();
	m_uiLength = 0;
}


/**********************************************************************
*	Biodata structure -- Conserved Domain information
***********************************************************************/
struct TOflCDInfo: public TOflClusterInfo
{
	double m_dMinBitScore;
	unsigned int m_uiHierarchyRoot;	//root pssmid
	unsigned int m_uiHierarchyParent;	//root pssmid
	unsigned int m_uiClusterPSSMID;
	bool m_bCurated;
	bool m_bIsStructDom;
	bool m_bMultiDom;
	
	list<TOflCDFeat> m_lstSpecFeatures;
	list<TOflCDFeat> m_lstGenFeatures;
	TOflCDInfo(void);
	void Reset(void);
};

TOflCDInfo::TOflCDInfo(void):
	TOflClusterInfo(), m_dMinBitScore(0.0),
	m_uiHierarchyRoot(0), m_uiHierarchyParent(0), m_uiClusterPSSMID(0), m_bCurated(false), m_bIsStructDom(false), m_bMultiDom(true), m_lstSpecFeatures(), m_lstGenFeatures()
{}


void TOflCDInfo::Reset(void)
{
	TOflClusterInfo::Reset();
	m_dMinBitScore = 0.0;
	m_uiHierarchyRoot = 0;
	m_uiHierarchyParent = 0;
	m_uiClusterPSSMID = 0;
	m_bCurated = false;
	m_bIsStructDom = false;
	m_bMultiDom = true;
	m_lstSpecFeatures.clear();
	m_lstGenFeatures.clear();
}

/**********************************************************************
*	Biodata structure -- Domain/Cluster info data center
***********************************************************************/
class COflDomClstInfo
{
public:
	
	static const unsigned int CDD_DATA_NOT_FOUND = 0x1;
	static const unsigned int CLUSTER_LINK_NOT_FOUND = 0x1 << 1;
	static const unsigned int HIERARCHY_DATA_NOT_FOUND = 0x1 << 2;
	static const unsigned int FEATURE_DATA_NOT_FOUND = 0x1 << 3;
	static const unsigned int SPTHRESHOLD_DATA_NOT_FOUND = 0x1 << 4;
	static const unsigned int GENERIC_FEATURE_DATA_NOT_FOUND = 0x1 << 5;
	COflDomClstInfo(void);
	//COflDomClstInfo(const CNcbiRegistry &reg);
	const TOflCDInfo* FindCDInfo(unsigned int pssmid) const;
	const TOflClusterInfo* FindClusterInfo(unsigned int pssmid) const;
	unsigned int GetStatus(void) const {return m_uiFlags;}
	void LoadData(const CNcbiRegistry &reg);	//load data and set m_uiFlags
	void LoadData(void);
private:
	struct TSortCDMapIterByAcc
	{
		bool operator () (map<unsigned int, TOflCDInfo> :: iterator v1, map<unsigned int, TOflCDInfo> :: iterator v2) {return v1->second.m_strAccession < v2->second.m_strAccession;}
	};
	
	static unsigned int GetPSSMIdByAcc(vector<map<unsigned int, TOflCDInfo> :: iterator> src, const string & acc);

	
	void x_LoadData(ifstream &cdid, ifstream &cdtrack, ifstream &clstlink, ifstream &cddannot, ifstream &cddannot_gen, ifstream &spthr);
	map<unsigned int, TOflCDInfo> m_mapOflCDInfo;
	map<unsigned int, TOflClusterInfo> m_mapOflClusterInfo;
		
	unsigned int m_uiFlags;
};

unsigned int COflDomClstInfo::GetPSSMIdByAcc(vector<map<unsigned int, TOflCDInfo> :: iterator> src, const string & acc)
{
	size_t rbegin = 0, rend = src.size(), dist = rend - rbegin;

	while (dist > 0)
	{
		size_t rmid = rbegin + (dist >> 1);

		if (src[rmid]->second.m_strAccession == acc) return src[rmid]->first;
		else if (src[rmid]->second.m_strAccession < acc) rbegin = rmid + 1;
		else rend = rmid;
		dist = rend - rbegin;
	}
	
	return 0;
}

COflDomClstInfo::COflDomClstInfo(void): m_mapOflCDInfo(), m_mapOflClusterInfo(), m_uiFlags(-1)
{
	//if (m_mapOflCDInfo.empty()) LoadData();
}

	
//COflDomClstInfo::COflDomClstInfo(const CNcbiRegistry &reg): m_mapOflCDInfo(), m_mapOflClusterInfo(), m_uiFlags(-1)
//{
//	//if (m_mapOflCDInfo.empty()) LoadData(reg);
//}

const TOflCDInfo* COflDomClstInfo::FindCDInfo(unsigned int pssmid) const
{
	map<unsigned int, TOflCDInfo> :: const_iterator iter = m_mapOflCDInfo.find(pssmid);
	if (m_mapOflCDInfo.end() == iter) return NULL;
	return &(iter->second);
}

const TOflClusterInfo* COflDomClstInfo::FindClusterInfo(unsigned int pssmid) const
{
	map<unsigned int, TOflClusterInfo> :: const_iterator iter = m_mapOflClusterInfo.find(pssmid);
	if (m_mapOflClusterInfo.end() == iter) return NULL;
	return &(iter->second);
}

void COflDomClstInfo::x_LoadData(ifstream &cdid, ifstream &cdtrack, ifstream &clstlink, ifstream &cddannot, ifstream &cddannot_gen, ifstream &spthr)
{
	m_mapOflCDInfo.clear();
	m_mapOflClusterInfo.clear();
	// -- first read cdid
	pair<unsigned int, TOflCDInfo> cdValue;
	pair<unsigned int, TOflClusterInfo> clstValue;
	
	vector<map<unsigned int, TOflCDInfo> :: iterator> vecAccSearchMap;
	vecAccSearchMap.reserve(81920);	//estimated.
	
	string buf;
	
	while (cdid.good())
	{
		getline(cdid, buf);
		if (!buf.empty() && '#' != buf[0])
		{

			size_t pos0 = 0, pos1 = buf.find('\t', pos0);
			unsigned int pssmid = (unsigned int)atoi(buf.substr(pos0, pos1 - pos0).c_str());
			
			pos0 = pos1 + 1;
			pos1 = buf.find('\t', pos0);
			
			
			if (('c' == buf[pos0] || 'C' == buf[pos0]) && ('l' == buf[pos0 + 1] || 'L' == buf[pos0 + 1]) && (buf[pos0 + 2] >= '0' && buf[pos0 + 2] <= '9'))	//cluster
			{
				clstValue.first = pssmid;
				map<unsigned int, TOflClusterInfo> :: iterator iterClst = m_mapOflClusterInfo.insert(clstValue).first;
				iterClst->second.m_strAccession = buf.substr(pos0, pos1 - pos0);
				
				pos0 = pos1 + 1;
				pos1 = buf.find('\t', pos0);
				iterClst->second.m_strShortName = buf.substr(pos0, pos1 - pos0);
				
				pos0 = pos1 + 1;
				pos1 = buf.find('\t', pos0);
				iterClst->second.m_strDefline = buf.substr(pos0, pos1 - pos0);
				
				pos0 = pos1 + 1;
				iterClst->second.m_uiLength = (unsigned int)atoi(buf.substr(pos0).c_str());
				continue;
			}
			
			cdValue.first = pssmid;
			map<unsigned int, TOflCDInfo> :: iterator iterDom = m_mapOflCDInfo.insert(cdValue).first;
			iterDom->second.m_strAccession = buf.substr(pos0, pos1 - pos0);
			
			if (iterDom->second.m_strAccession[2] >= '0' && iterDom->second.m_strAccession[2] <= '9')	//cdxxx, sdxxx or clxxx
			{
				if ('d' == iterDom->second.m_strAccession[1] || 'D' == iterDom->second.m_strAccession[1])	//cd or sd
				{
					if ('c' == iterDom->second.m_strAccession[0] || 'C' == iterDom->second.m_strAccession[0])	//curated
					{
						iterDom->second.m_bCurated = true;
						vecAccSearchMap.push_back(iterDom);
					}
						
					else if ('s' == iterDom->second.m_strAccession[0] || 'S' == iterDom->second.m_strAccession[0])	//curated
						iterDom->second.m_bIsStructDom = true;
				}
			}
			
			
			pos0 = pos1 + 1;
			pos1 = buf.find('\t', pos0);
			iterDom->second.m_strShortName = buf.substr(pos0, pos1 - pos0);
			
			pos0 = pos1 + 1;
			pos1 = buf.find('\t', pos0);
			iterDom->second.m_strDefline = buf.substr(pos0, pos1 - pos0);
			
			pos0 = pos1 + 1;
			iterDom->second.m_uiLength = (unsigned int)atoi(buf.substr(pos0).c_str());
		}
	}
	
	
	sort(vecAccSearchMap.begin(), vecAccSearchMap.end(), TSortCDMapIterByAcc());
	map<unsigned int, TOflCDInfo> :: iterator iterMapEnd = m_mapOflCDInfo.end();
	
	while (cdtrack.good())
	{
		getline(cdtrack, buf);
		if (!buf.empty() && '#' != buf[0])
		{
			// -- accession
			size_t pos0 = 0, pos1 = buf.find(' ', pos0);
			string currAcc = buf.substr(pos0, pos1 - pos0);
			
			// -- title
			while(' ' == buf[pos1]) ++pos1;
			pos0 = pos1;
			pos1 = buf.find(' ', pos0);
			
			
			// -- pssmid
			while(' ' == buf[pos1]) ++pos1;
			pos0 = pos1;
			pos1 = buf.find(' ', pos0);
			
			unsigned int pssmid = (unsigned int)atoi(buf.substr(pos0, pos1 - pos0).c_str());
			
			map<unsigned int, TOflCDInfo> :: iterator iterDom = m_mapOflCDInfo.find(pssmid);
			if (iterMapEnd != iterDom)
			{

				// -- parent acc
				while(' ' == buf[pos1]) ++pos1;
				pos0 = pos1;
				pos1 = buf.find(' ', pos0);
				
				string parentAcc = buf.substr(pos0, pos1 - pos0);
				if ("N/A" != parentAcc)
				{
					if (parentAcc == iterDom->second.m_strAccession)
						iterDom->second.m_uiHierarchyParent = iterDom->first;
					else
						iterDom->second.m_uiHierarchyParent = GetPSSMIdByAcc(vecAccSearchMap, parentAcc);
				}
				// -- root acc
				while(' ' == buf[pos1]) ++pos1;
				pos0 = pos1;
				pos1 = buf.find(' ', pos0);
				string rootAcc = buf.substr(pos0, pos1 - pos0);
				
				if (rootAcc == parentAcc)
					iterDom->second.m_uiHierarchyRoot = iterDom->second.m_uiHierarchyParent;
				else if (rootAcc == iterDom->second.m_strAccession)
					iterDom->second.m_uiHierarchyRoot = iterDom->first;
				else
					iterDom->second.m_uiHierarchyRoot = GetPSSMIdByAcc(vecAccSearchMap, rootAcc);
			}
		}
	}	//cdtrack
	
	

	while (clstlink.good())
	{
		getline(clstlink, buf);
		if (!buf.empty() && '#' != buf[0])
		{
			size_t pos0 = 0, pos1 = buf.find('\t', pos0);
			pos0 = pos1 + 1;
			pos1 = buf.find('\t', pos0);
			
			unsigned int pssmid = (unsigned int)atoi(buf.substr(pos0, pos1 - pos0).c_str());
			map<unsigned int, TOflCDInfo> :: iterator iterDom = m_mapOflCDInfo.find(pssmid);
  
			if (iterMapEnd != iterDom)
			{
				pos0 = pos1 + 1;
				pos1 = buf.find('\t', pos0);
				pos0 = pos1 + 1;
				//pos1 = buf.find('\t', pos0);
				
				iterDom->second.m_uiClusterPSSMID = (unsigned int)atoi(buf.substr(pos0).c_str());
				iterDom->second.m_bMultiDom = false;
			}
		}
	}
	
	
	
	TOflCDFeat featValue;


	while (cddannot.good())
	{
		getline(cddannot, buf);
		if (!buf.empty() && '#' != buf[0])
		{

			size_t pos0 = 0, pos1 = buf.find('\t', pos0);
			
			unsigned int pssmid = (unsigned int)atoi(buf.substr(pos0, pos1 - pos0).c_str());
			
			
			map<unsigned int, TOflCDInfo> :: iterator iterDom = m_mapOflCDInfo.find(pssmid);


			
			if (iterMapEnd != iterDom)
			{
				iterDom->second.m_lstSpecFeatures.push_back(featValue);
				TOflCDFeat &tgt =  *(iterDom->second.m_lstSpecFeatures.rbegin());
				
				pos0 = pos1 + 1;
				pos1 = buf.find('\t', pos0);	//accession
				
				pos0 = pos1 + 1;
				pos1 = buf.find('\t', pos0);	//shortname
				
				pos0 = pos1 + 1;
				pos1 = buf.find('\t', pos0);	//index
				
				tgt.m_iIndex = atoi(buf.substr(pos0, pos1 - pos0).c_str());
				
				pos0 = pos1 + 1;
				pos1 = buf.find('\t', pos0);	//title
				
				tgt.m_strTitle = buf.substr(pos0, pos1 - pos0);
				
				pos0 = pos1 + 1;
				pos1 = buf.find('\t', pos0);	//Motif
				
				if ('0' != buf[pos0])
					tgt.SetMotifStr(buf.substr(pos0, pos1 - pos0).c_str());
				
				pos0 = pos1 + 1;
				pos1 = buf.find('\t', pos0);	//structure based evidence present
				
				if ('1' == buf[pos0])
					tgt.flags |= TOflCDFeat::STRUCTURE_BASED_EVIDENCE;
					
				pos0 = pos1 + 1;
				pos1 = buf.find('\t', pos0);	//reference based evidence present
				
				if ('1' == buf[pos0])
					tgt.flags |= TOflCDFeat::REFERENCE_BASED_EVIDENCE;
					
				pos0 = pos1 + 1;
				pos1 = buf.find('\t', pos0);	//additional comments present
				
				if ('1' == buf[pos0])
					tgt.flags |= TOflCDFeat::ADDITIONAL_COMMENTS;
					
				pos0 = pos1 + 1;
				pos1 = buf.find('\t', pos0);	//site coordinates
				
				vector<int> vecSitePos;	//to store coordinates
				
				vecSitePos.reserve((pos1 - pos0) * 3);	//estimated
				
				size_t pos_c = buf.find(',', pos0);
				
				while (pos_c < pos1)
				{
					vecSitePos.push_back(atoi(buf.substr(pos0, pos_c - pos0).c_str()));
					pos0 = pos_c + 1;
					pos_c = buf.find(',', pos0);
				}
				
				// -- last one
				vecSitePos.push_back(atoi(buf.substr(pos0, pos1 - pos0).c_str()));
				
				tgt.SetData(vecSitePos);
				
				pos0 = pos1 + 1;
				pos1 = buf.find('\t', pos0);	//type
				
				if (string::npos == pos1)	//no other field
					tgt.m_iType = atoi(buf.substr(pos0).c_str());
				else	//assume has motif string added
					tgt.m_iType = atoi(buf.substr(pos0, pos1 - pos0).c_str());
			}
		}
	}

	while (cddannot_gen.good())
	{
		getline(cddannot_gen, buf);
		if (!buf.empty() && '#' != buf[0])
		{

			size_t pos0 = 0, pos1 = buf.find('\t', pos0);
			
			unsigned int pssmid = (unsigned int)atoi(buf.substr(pos0, pos1 - pos0).c_str());
			
			
			map<unsigned int, TOflCDInfo> :: iterator iterDom = m_mapOflCDInfo.find(pssmid);
			
			if (iterMapEnd != iterDom)
			{
				iterDom->second.m_lstGenFeatures.push_back(featValue);
				TOflCDFeat &tgt =  *(iterDom->second.m_lstGenFeatures.rbegin());
				
				pos0 = pos1 + 1;
				pos1 = buf.find('\t', pos0);	//accession
				
				pos0 = pos1 + 1;
				pos1 = buf.find('\t', pos0);	//shortname
				
				pos0 = pos1 + 1;
				pos1 = buf.find('\t', pos0);	//index
				
				tgt.m_iIndex = atoi(buf.substr(pos0, pos1 - pos0).c_str());
				
				pos0 = pos1 + 1;
				pos1 = buf.find('\t', pos0);	//title
				
				tgt.m_strTitle = buf.substr(pos0, pos1 - pos0);
				
				pos0 = pos1 + 1;
				pos1 = buf.find('\t', pos0);	//Motif
				
				if ('0' != buf[pos0])
					tgt.SetMotifStr(buf.substr(pos0, pos1 - pos0).c_str());
				
				pos0 = pos1 + 1;
				pos1 = buf.find('\t', pos0);	//structure based evidence present
				
				if ('1' == buf[pos0])
					tgt.flags |= TOflCDFeat::STRUCTURE_BASED_EVIDENCE;
					
				pos0 = pos1 + 1;
				pos1 = buf.find('\t', pos0);	//reference based evidence present
				
				if ('1' == buf[pos0])
					tgt.flags |= TOflCDFeat::REFERENCE_BASED_EVIDENCE;
					
				pos0 = pos1 + 1;
				pos1 = buf.find('\t', pos0);	//additional comments present
				
				if ('1' == buf[pos0])
					tgt.flags |= TOflCDFeat::ADDITIONAL_COMMENTS;
					
				pos0 = pos1 + 1;
				pos1 = buf.find('\t', pos0);	//site coordinates
				
				vector<int> vecSitePos;	//to store coordinates
				
				vecSitePos.reserve((pos1 - pos0) * 3);	//estimated
				
				size_t pos_c = buf.find(',', pos0);
				
				while (pos_c < pos1)
				{
					vecSitePos.push_back(atoi(buf.substr(pos0, pos_c - pos0).c_str()));
					pos0 = pos_c + 1;
					pos_c = buf.find(',', pos0);
				}
				
				// -- last one
				vecSitePos.push_back(atoi(buf.substr(pos0, pos1 - pos0).c_str()));
				
				tgt.SetData(vecSitePos);
				
				pos0 = pos1 + 1;
				pos1 = buf.find('\t', pos0);	//type
				
				if (string::npos == pos1)	//no other field
					tgt.m_iType = atoi(buf.substr(pos0).c_str());
				else	//assume has motif string added
					tgt.m_iType = atoi(buf.substr(pos0, pos1 - pos0).c_str());
			}
		}
	}
		
	

	while (spthr.good())
	{
		getline(spthr, buf);

		if (!buf.empty() && '#' != buf[0])
		{

			size_t pos0 = 0, pos1 = buf.find('\t', pos0);

			unsigned int pssmid = (unsigned int)atoi(buf.substr(pos0, pos1 - pos0).c_str());

			
			map<unsigned int, TOflCDInfo> :: iterator iterDom = m_mapOflCDInfo.find(pssmid);
			
			if (iterMapEnd != iterDom)
			{
				// -- accession, no need
				pos0 = pos1 + 1;
				pos1 = buf.find('\t', pos0);

				
				// -- score
				pos0 = pos1 + 1;
				pos1 = buf.find('\t', pos0);


				if (string::npos == pos1)
				{

					iterDom->second.m_dMinBitScore = atof(buf.substr(pos0).c_str());
				}
				else
				{


					iterDom->second.m_dMinBitScore = atof(buf.substr(pos0, pos1 - pos0).c_str());
				}
			}
		}
	}
}
void COflDomClstInfo::LoadData(void)
{
	m_uiFlags = 0;
	ifstream cdid, cdtrack, clstlink, cddannot, cddannot_gen, spthr;
	
	if (!g_strDataPath.empty())
	{
		cdid.open((g_strDataPath + "/" + CDIDFILE).c_str());
		if (!cdid.good())
			m_uiFlags |= CDD_DATA_NOT_FOUND;
			
		cdtrack.open((g_strDataPath + "/" + CDTRACKFILE).c_str());
		if (!cdtrack.good())
			m_uiFlags |= HIERARCHY_DATA_NOT_FOUND;
			
		clstlink.open((g_strDataPath + "/" + CLSTLINKFILE).c_str());
		if (!clstlink.good())
			m_uiFlags |= CLUSTER_LINK_NOT_FOUND;
			
		cddannot.open((g_strDataPath + "/" + SPFEATFILE).c_str());
		if (!cddannot.good())
			m_uiFlags |= FEATURE_DATA_NOT_FOUND;
		
		cddannot_gen.open((g_strDataPath + "/" + GENFEATFILE).c_str());
		if (!cddannot_gen.good())
			m_uiFlags |= GENERIC_FEATURE_DATA_NOT_FOUND;
			
		spthr.open((g_strDataPath + "/" + MINBSCOREFILE).c_str());
		if (!spthr.good())
			m_uiFlags |= SPTHRESHOLD_DATA_NOT_FOUND;
	}
	else
	{
		cdid.open(CDIDFILE);
		if (!cdid.good())
			m_uiFlags |= CDD_DATA_NOT_FOUND;
			
		cdtrack.open(CDTRACKFILE);
		if (!cdtrack.good())
			m_uiFlags |= HIERARCHY_DATA_NOT_FOUND;
			
		clstlink.open(CLSTLINKFILE);
		if (!clstlink.good())
			m_uiFlags |= CLUSTER_LINK_NOT_FOUND;
			
		cddannot.open(SPFEATFILE);
		if (!cddannot.good())
			m_uiFlags |= FEATURE_DATA_NOT_FOUND;
		
		cddannot_gen.open(GENFEATFILE);
		if (!cddannot_gen.good())
			m_uiFlags |= GENERIC_FEATURE_DATA_NOT_FOUND;
			
		spthr.open(MINBSCOREFILE);
		if (!spthr.good())
			m_uiFlags |= SPTHRESHOLD_DATA_NOT_FOUND;
	}
	
	if (m_uiFlags > 0)
	{
		cdid.close();
		cdtrack.close();
		clstlink.close();
		cddannot.close();
		cddannot_gen.close();
		spthr.close();
		return;
	}
	
	x_LoadData(cdid, cdtrack, clstlink, cddannot, cddannot_gen, spthr);
	
	cdid.close();
	cdtrack.close();
	clstlink.close();
	cddannot.close();
	cddannot_gen.close();
	spthr.close();
}

void COflDomClstInfo::LoadData(const CNcbiRegistry &reg)
{

	// -- check everything and see if ok                                
	m_uiFlags = 0;
	
	ifstream cdid, cdtrack, clstlink, cddannot, cddannot_gen, spthr;
	
	if (!g_strDataPath.empty())
	{
		cdid.open((g_strDataPath + "/" + CDIDFILE).c_str());
		if (!cdid.good())
			m_uiFlags |= CDD_DATA_NOT_FOUND;
			
		cdtrack.open((g_strDataPath + "/" + CDTRACKFILE).c_str());
		if (!cdtrack.good())
			m_uiFlags |= HIERARCHY_DATA_NOT_FOUND;
			
		clstlink.open((g_strDataPath + "/" + CLSTLINKFILE).c_str());
		if (!clstlink.good())
			m_uiFlags |= CLUSTER_LINK_NOT_FOUND;
			
		cddannot.open((g_strDataPath + "/" + SPFEATFILE).c_str());
		if (!cddannot.good())
			m_uiFlags |= FEATURE_DATA_NOT_FOUND;
		
		cddannot_gen.open((g_strDataPath + "/" + GENFEATFILE).c_str());
		if (!cddannot_gen.good())
			m_uiFlags |= GENERIC_FEATURE_DATA_NOT_FOUND;
		
		spthr.open((g_strDataPath + "/" + MINBSCOREFILE).c_str());
		if (!spthr.good())
			m_uiFlags |= SPTHRESHOLD_DATA_NOT_FOUND;
	}
	else
	{
		cdid.open(reg.Get(DATAPATH, CDDIDS).c_str());
		if (!cdid.good())
			m_uiFlags |= CDD_DATA_NOT_FOUND;
			
		cdtrack.open(reg.Get(DATAPATH, CDTRACKINFO).c_str());
		if (!cdtrack.good())
			m_uiFlags |= HIERARCHY_DATA_NOT_FOUND;
			
		clstlink.open(reg.Get(DATAPATH, CLUSTERLINKS).c_str());
		if (!clstlink.good())
			m_uiFlags |= CLUSTER_LINK_NOT_FOUND;
		
		cddannot.open(reg.Get(DATAPATH, FEATURES).c_str());
		if (!cddannot.good())
			m_uiFlags |= FEATURE_DATA_NOT_FOUND;
		
		cddannot_gen.open(reg.Get(DATAPATH, GENERIC_FEATURES).c_str());
		if (!cddannot_gen.good())
			m_uiFlags |= GENERIC_FEATURE_DATA_NOT_FOUND;
		
		spthr.open(reg.Get(DATAPATH, SPECIFICTHRESHOLDS).c_str());
		if (!spthr.good())
			m_uiFlags |= SPTHRESHOLD_DATA_NOT_FOUND;
	}
	
	if (m_uiFlags > 0)
	{
		cdid.close();
		cdtrack.close();
		clstlink.close();
		cddannot.close();
		cddannot_gen.close();
		spthr.close();
		return;
	}
	
	x_LoadData(cdid, cdtrack, clstlink, cddannot, cddannot_gen, spthr);
	
	cdid.close();
	cdtrack.close();
	clstlink.close();
	cddannot.close();
	cddannot_gen.close();
	spthr.close();
	
}


/**********************************************************************
*	Biodata structure -- sequence alignments
***********************************************************************/
struct TOflAlignInfo
{
	enum EIndex
	{
		eEnumStart = 0,
		// ------------------
		SORT_BY_EVALUE = eEnumStart,
		SORT_BY_BITSCORE,
		SORT_BY_SEQ_IDENTITY,
		SORT_BY_ALIGNED_LENGTH,
		
		// ------------------
		eEnumStop
	};
	
	
	static const EIndex eDefault = SORT_BY_EVALUE;
	static const char* dimLits[eEnumStop];
	static const char* dimLabels[eEnumStop];
	
	typedef bool LPFN_COMPARE(const TOflAlignInfo *p1, const TOflAlignInfo *p2);
	
	static bool EValueCompare(const TOflAlignInfo *p1, const TOflAlignInfo *p2)
	{
		return p1->m_dEValue < p2->m_dEValue;
	}
	
	static bool BitScoreCompare(const TOflAlignInfo *p1, const TOflAlignInfo *p2)
	{
		return p1->m_dBitScore > p2->m_dBitScore;
	}
	
	static bool AlignedLengthCompare(const TOflAlignInfo *p1, const TOflAlignInfo *p2)
	{
		return p1->m_uiAlignedLen > p2->m_uiAlignedLen;
	}
	
	static bool SeqIdentityCompare(const TOflAlignInfo *p1, const TOflAlignInfo *p2)
	{
		return p1->m_dSeqIdentity > p2->m_dSeqIdentity;
	}
	
	struct TSortObj
	{
		TSortObj(LPFN_COMPARE lpfnCompare): m_lpfnCompare(lpfnCompare) {};
		TSortObj(int iSortIdx);
		LPFN_COMPARE * GetCompareFunc(void) const {return m_lpfnCompare;}
		
		bool operator () (const TOflAlignInfo *p1, const TOflAlignInfo *p2)
		{
			return m_lpfnCompare(p1, p2);
		}
		
		LPFN_COMPARE *m_lpfnCompare;
	};
	
	// -- scores
	unsigned int m_uiAlignedLen;
	double m_dAlignedPct;
	int m_iScore;
	double m_dEValue;
	double m_dBitScore;
	int m_iNumIdent;
	double m_dSeqIdentity;
	
	// -- added 6/27/2011 -- for na aligns
	enum EAlignType
	{
		eNormal,	//Na to Na or prot to prot
		ePr2Na
	} m_eAlignType;
	// -- added 6-23-2011: for NA search
	//enum ENa_strand NCBI_PACKED_ENUM_TYPE( unsigned char )
	//{
    //	eNa_strand_unknown  =   0,
    //	eNa_strand_plus     =   1,
    //	eNa_strand_minus    =   2,
    //	eNa_strand_both     =   3,  ///< in forward orientation
    //	eNa_strand_both_rev =   4,  ///< in reverse orientation
    //	eNa_strand_other    = 255
	//} NCBI_PACKED_ENUM_END();
	ENa_strand m_eStrand;
	int m_iReadingFrame;	//seqlen << 2 | (abs(rf) - 1)
	//int m_iSlaveCoordConvert;	//a slave coordination conversion factor for negative factor
	
	// -- processed align info
	vector<int> m_vecMStarts;
	vector<int> m_vecSStarts;
	vector<unsigned int> m_vecLens;
	int m_iRegionIdx;	//region on the query sequence. 
	
	CSegSet m_ClipSet;
	// -- methods
	TOflAlignInfo(void);
	
	
	void Pr2NaConvert(CSegSet &segset) const;
	void MapSegSet(CSegSet &segset, bool doConvert = true) const;
	void CalcMasterGaps(unsigned int gapThreshold, CSegSet &segset) const;
	
	int GetRFIdx(void) const;
	
	int NAPlus2Pr(int na) const {return na / RF_SIZE;}
	int Pr2NAPlus(int pr) const	{return pr * RF_SIZE + (m_iReadingFrame & RF_SIZE);}
	int NAMinus2Pr(int na) const {return ((m_iReadingFrame >> 2) - na - 1) / RF_SIZE;}
	int Pr2NAMinus(int pr) const {return (m_iReadingFrame >> 2) - pr * RF_SIZE - (m_iReadingFrame & RF_SIZE) - 1;}
	
	//protein coord in pr, return true if plus strand, false if minus strand
	bool NA2Pr(int na, int &pr) const;
	bool Pr2NA(int pr, int &na) const;

};

const char * TOflAlignInfo::dimLits[] = {"evalue", "bitscore", "seqidentity", "alignedlen"};
const char * TOflAlignInfo::dimLabels[] = {"BLAST E-value", "BLAST bit score", "Sequence Identity", "Aligned Length"};
TOflAlignInfo::TSortObj::TSortObj(int iSortIdx):
	m_lpfnCompare(TOflAlignInfo::EValueCompare)
{
	switch (iSortIdx)
	{
		case TOflAlignInfo::SORT_BY_BITSCORE:
			m_lpfnCompare = TOflAlignInfo::BitScoreCompare;
			break;
		case TOflAlignInfo::SORT_BY_SEQ_IDENTITY:
			m_lpfnCompare = TOflAlignInfo::SeqIdentityCompare;
			break;
		case TOflAlignInfo::SORT_BY_ALIGNED_LENGTH:
			m_lpfnCompare = TOflAlignInfo::AlignedLengthCompare;
			break;
		default:;
	}
}

TOflAlignInfo::TOflAlignInfo(void):
	m_uiAlignedLen(0), m_dAlignedPct(0.0), m_iScore(0), m_dEValue(0.0), m_dBitScore(0.0), m_iNumIdent(0), m_dSeqIdentity(0.0), m_eAlignType(eNormal), m_eStrand(eNa_strand_unknown), m_iReadingFrame(0), /*m_iSlaveCoordConvert(-1), */m_vecMStarts(), m_vecSStarts(), m_vecLens(), m_iRegionIdx(0), m_ClipSet()
{}

void TOflAlignInfo::Pr2NaConvert(CSegSet &segset) const
{

	if (ePr2Na == m_eAlignType)
	{
		
		segset.m_iFactor *= RF_SIZE;
		
		if ((m_iReadingFrame >> 2) > 0)	//is minus
		{
			for (TSegs::iterator iterSeg = segset.m_lstContainer.begin(); iterSeg != segset.m_lstContainer.end(); ++iterSeg)
			{
				
				int newFrom = Pr2NAMinus(iterSeg->from);
				iterSeg->from = Pr2NAMinus(iterSeg->to) - RF_SIZE + 1;
				iterSeg->to = newFrom;
				EResidueExtension oriLext = iterSeg->lext;
				iterSeg->lext = iterSeg->rext;
				iterSeg->rext = oriLext;
				
				// -- reverse ori_as well
				iterSeg->ori_from += (iterSeg->to - iterSeg->from) / segset.m_iFactor;
				
			}
			segset.m_iFactor = -segset.m_iFactor;
			segset.m_lstContainer.reverse();
		}
		else	//plus strand
		{
			for (TSegs::iterator iter = segset.m_lstContainer.begin(); iter != segset.m_lstContainer.end(); ++iter)
			{
				iter->from = Pr2NAPlus(iter->from);
				iter->to = Pr2NAPlus(iter->to) + RF_SIZE - 1;
			}
		}
	}
}

void TOflAlignInfo::MapSegSet(CSegSet &segset, bool doConvert) const
{
	if (segset.IsEmpty()) return;
	
	TSegs::iterator iter = segset.m_lstContainer.begin();
	size_t idx = 0;
	
	while (iter != segset.m_lstContainer.end())
	{
		if (idx >= m_vecLens.size() || m_vecSStarts[idx] > iter->to)	// discard this seg
		{
			TSegs::iterator temp = iter;
			++iter;
			segset.m_lstContainer.erase(temp);
		}
		else
		{
			int diff = m_vecSStarts[idx] - m_vecMStarts[idx];
			int end = m_vecSStarts[idx] + m_vecLens[idx] - 1;
			if (end >= iter->to)
			{
				if (iter->from < m_vecSStarts[idx])	//segment shrinked from left. deal with ori_from
				{
					iter->ori_from += (m_vecSStarts[idx] - iter->from) / segset.m_iFactor;
					iter->from = m_vecSStarts[idx];
				}
			
				if (end == iter->to)	//seg happens to end at aligned seg
					++idx;	//here is the chance to advance idx
				
				
				// -- mapping
				iter->from -= diff;
				iter->to -= diff;
				
				++iter;
			}
			else if (end >= iter->from)	//end < iter->to
			{
				TSeg temp(iter->from, end, iter->lext, eAccurate);
				temp.ori_from = iter->ori_from;
				
				if (temp.from < m_vecSStarts[idx])
				{
					temp.ori_from += (m_vecSStarts[idx] - temp.from) / segset.m_iFactor;
					temp.from = m_vecSStarts[idx];
					temp.lext = eAccurate;
				}
				
				// -- mapping
				temp.from -= diff;
				temp.to -= diff;
				segset.m_lstContainer.insert(iter, temp);
				
				// -- cut original seg for next round
				iter->ori_from += (end - iter->from + 1) / segset.m_iFactor;
				iter->from = end + 1;
				//iter->ori_from += (iter->from) / segset.m_iFactor;
				iter->lext = eAccurate;
				++idx;
			}
			else	//end < iter->from, step to next denseg
				++idx;
		}
	}
	if (doConvert) Pr2NaConvert(segset);
}

void TOflAlignInfo::CalcMasterGaps(unsigned int gapThreshold, CSegSet &segset) const
{
	segset.Clear();
	if (!m_vecLens.empty())
	{
		size_t segs = m_vecLens.size();
		for (size_t i = 0; i < segs - 1; ++i)
		{
			int gapstart = m_vecMStarts[i] + (int)m_vecLens[i];
			int gaplen = m_vecMStarts[i + 1] - gapstart;
			
			if (gaplen >= (int)gapThreshold)	//consider a gap
				segset.AddSeg(gapstart, gapstart + gaplen - 1);
		}
	}
}

int TOflAlignInfo::GetRFIdx(void) const
{
	int rfAbs = m_iReadingFrame & RF_SIZE;	//0, 1, 2
	if (m_iReadingFrame >> 2 > 0)	//minus chain
	{

		return rfAbs + RF_SIZE;
	}
	return rfAbs;
}

bool TOflAlignInfo::NA2Pr(int na, int &pr) const
{
	int len = m_iReadingFrame >> 2;
	if (len > 0)	//minus strand
	{
		pr = (len - na - 1) / 3;
		return false;
	}
	pr = na / 3;
	return true;
}

bool TOflAlignInfo::Pr2NA(int pr, int &na) const
{
	int len = m_iReadingFrame >> 2;
	if (len > 0)
	{
		na = len - pr * 3 - (m_iReadingFrame & RF_SIZE) - 1;
		return false;
	}
	na = pr * 3 + (m_iReadingFrame & RF_SIZE);
	return true;
}



struct TOflCDAlignInfo: public TOflAlignInfo
{
	static const TSeqPos GAP_THRESHOLD = 35;
	unsigned int m_uiPSSMID;
	// -- calculated
	int m_iFrom;
	int m_iTo;
	double m_dNMissing;
	double m_dCMissing;
	bool m_bSpecQualified;	//higher bitscore than threshold
	int m_iRepClass;	//single and multi -- sort separately
	bool m_bRep;
	bool m_bLifted;	//lifted by a higher evalue and approved by architecture frequencs
	
	TOflCDAlignInfo(void);
	
	//virtual CJSVar CreateJsonNode(int opts = 0) const;
	
	// -- including properly mapped features
	
	// -- extra
	bool IsSpecific(void) const {return m_bRep && m_bSpecQualified;}
	
	// -- map a feature to dst. return 0: success. Return 1: not all
	// -- essential residues are mapped. return 2: residue type
	// -- mismatch (if seqData is not empty)
	
	// -- convert segs (master coordinates already) to translated cooridinates so TOflCDFeat can perform motif check.
	// -- return reading frame for the alignment, which can be used to select a translation frame index (0 - 5). if 
	// -- return 0 for protein sequence that needs no translation.
	int GetTranslatedPosMap(const CSegSet &segs, vector<CSegSet::TResiduePos> &rMappedAAPos) const;
};

TOflCDAlignInfo::TOflCDAlignInfo(void):
	TOflAlignInfo(), m_uiPSSMID(0), m_iFrom(0), m_iTo(0), m_dNMissing(0.0), m_dCMissing(0.0), m_bSpecQualified(false), m_iRepClass(0), m_bRep(false), m_bLifted(false)
{}



// -- segs are already mapped to master sequence. this is to convert to traslated coordinates. return reading frame index (0-5)
int TOflCDAlignInfo::GetTranslatedPosMap(const CSegSet &segs, vector<CSegSet::TResiduePos> &rMappedAAPos) const
{
	rMappedAAPos.clear();
	segs.GetTranslatedPosMap(m_iReadingFrame >> 2, rMappedAAPos);
	if (ePr2Na == m_eAlignType)	//needs translation
		return GetRFIdx();
	else
		return 0;
	
	
	//rMappedAAPos.clear();
	//if (ePr2Na == m_eAlignType)	//needs translation
	//{
	//	segs.GetOriSimplePosMap(rMappedAAPos);
	//	int seqLen = m_iReadingFrame >> 2;
	//	if (seqLen > 0)
	//	{
	//		for (size_t i = 0, len = rMappedAAPos.size(); i < len; ++i)
	//			rMappedAAPos[i].curr = (seqLen - rMappedAAPos[i].curr - 1) / RF_SIZE;
	//		return GetRFIdx();	//isplus = false
	//	}
	//	else
	//	{
	//		for (size_t i = 0, len = rMappedAAPos.size(); i < len; ++i)
	//			rMappedAAPos[i].curr /= RF_SIZE;
	//		return GetRFIdx();	//isplus = true
	//	}
	//}
	//else
	//	segs.GetSimplePosMap(rMappedAAPos);
	//	
	//return 0;
		
}

/**********************************************************************
*	Biodata structure -- Align index
***********************************************************************/
struct TOflAlignIndice
{
	struct __TOflAlignRecord
	{
		const TOflCDAlignInfo * pAlign;
		const TOflCDInfo * pCDInfo;
		const TOflClusterInfo * pClst;
		const TOflCDInfo * pRootCDInfo;
		
		//const CDomainInfoMaps::TCuratedClusterInfo * pCuratedClst;
		size_t idx;	//need to give area id an index to align info
		size_t idxidx;	//the index in concise
		
		__TOflAlignRecord(void): pAlign(nullptr), pCDInfo(nullptr), pClst(nullptr), pRootCDInfo(nullptr), idx(-1), idxidx(-1) {};
	};
	
	
	vector<size_t> m_vecSortedIndice;
	vector<size_t> m_vecConciseIndice;
	vector<size_t> m_vecStdIndice;
	// -- modified 5/8/2012 -- now feature and rep hits are separated -- since we introduced non-NCBI-curated specific hits
	vector<size_t> m_vecQualifiedFeatIndice;	//should be every region's best evalue curated
	// -- modified 9/8/2014 -- Structure domains to add motif annotations
	vector<size_t> m_vecSDIndice;
	
	TOflAlignIndice(void): m_vecSortedIndice(), m_vecConciseIndice(), m_vecStdIndice(), m_vecQualifiedFeatIndice(), m_vecSDIndice() {};

	void CreateRecordSets(const vector<TOflCDAlignInfo> &rAlignments, const COflDomClstInfo & rDomInfo, vector<__TOflAlignRecord> &rDomAligns, vector<__TOflAlignRecord> &rFeatAligns, int mode) const;
	void ExtractFeatAligns(const __TOflAlignRecord &rRepRec, const vector<TOflCDAlignInfo> &rAlignments, const COflDomClstInfo & rDomInfo, vector<__TOflAlignRecord> &rResult) const;
};

void TOflAlignIndice::CreateRecordSets(const vector<TOflCDAlignInfo> &rAlignments, const COflDomClstInfo & rDomInfo, vector<__TOflAlignRecord> &rDomAligns, vector<__TOflAlignRecord> &rFeatAligns, int mode) const
{
	rDomAligns.clear();
	rFeatAligns.clear();
	
	size_t featIdx0 = 0;
	size_t amendCount = 0;
	size_t ccsBase = m_vecConciseIndice.size();


	if (TDataModes::e_rep == mode)	//all rep
	{
		for (size_t iidx = 0, iidxEnd = m_vecConciseIndice.size(); iidx < iidxEnd; ++iidx)
		{
			__TOflAlignRecord alignRec;
			alignRec.pAlign = &(rAlignments[m_vecConciseIndice[iidx]]);
			alignRec.pCDInfo = rDomInfo.FindCDInfo(alignRec.pAlign->m_uiPSSMID);
			if (alignRec.pCDInfo->m_uiClusterPSSMID > 0) alignRec.pClst = rDomInfo.FindClusterInfo(alignRec.pCDInfo->m_uiClusterPSSMID);
			if (alignRec.pCDInfo->m_uiHierarchyRoot > 0)
			{
				if (alignRec.pCDInfo->m_uiHierarchyRoot == alignRec.pAlign->m_uiPSSMID) alignRec.pRootCDInfo = alignRec.pCDInfo;
				else alignRec.pRootCDInfo = rDomInfo.FindCDInfo(alignRec.pCDInfo->m_uiHierarchyRoot);
			}
			alignRec.idx = m_vecConciseIndice[iidx];
			alignRec.idxidx = iidx;

			if (0 == alignRec.pAlign->m_iRepClass)	//non-multi
			{
				ExtractFeatAligns(alignRec, rAlignments, rDomInfo, rFeatAligns);
				
				size_t featIdx1 = rFeatAligns.size();
		    	
				for (size_t i = featIdx0; i < featIdx1; ++i)
				{
					if (rFeatAligns[i].idx == alignRec.idx)
						rFeatAligns[i].idxidx = alignRec.idxidx;
					else
						rFeatAligns[i].idxidx = iidxEnd + amendCount++;
				}
				
				featIdx0 = featIdx1;
			}
			
			rDomAligns.push_back(alignRec);
		}
	}
	else	//non-concise 
	{
		size_t repIdx = 0;
		
		size_t featIdx0 = 0;
		size_t amendCount = 0;
		
		const vector<size_t> &nonConciseIdx = (TDataModes::e_full == mode ? m_vecSortedIndice : m_vecStdIndice);
		
		for (size_t iidx = 0, iidxEnd = nonConciseIdx.size(); iidx < iidxEnd; ++iidx)
		{
			__TOflAlignRecord alignRec;
			alignRec.pAlign = &(rAlignments[nonConciseIdx[iidx]]);
			alignRec.pCDInfo = rDomInfo.FindCDInfo(alignRec.pAlign->m_uiPSSMID);


			if (alignRec.pCDInfo->m_uiClusterPSSMID > 0) alignRec.pClst = rDomInfo.FindClusterInfo(alignRec.pCDInfo->m_uiClusterPSSMID);
			if (alignRec.pCDInfo->m_uiHierarchyRoot > 0)
			{
				if (alignRec.pCDInfo->m_uiHierarchyRoot == alignRec.pAlign->m_uiPSSMID) alignRec.pRootCDInfo = alignRec.pCDInfo;
				else alignRec.pRootCDInfo = rDomInfo.FindCDInfo(alignRec.pCDInfo->m_uiHierarchyRoot);
			}
			alignRec.idx = nonConciseIdx[iidx];
			alignRec.idxidx = -1;
			if (alignRec.pAlign->m_bRep)
			{
				alignRec.idxidx = repIdx++;
				if (0 == alignRec.pAlign->m_iRepClass)	//non-multi
				{
					ExtractFeatAligns(alignRec, rAlignments, rDomInfo, rFeatAligns);
					size_t featIdx1 = rFeatAligns.size();
		    	
					for (size_t i = featIdx0; i < featIdx1; ++i)
					{
						if (rFeatAligns[i].idx == alignRec.idx)
							rFeatAligns[i].idxidx = alignRec.idxidx;
						else
							rFeatAligns[i].idxidx = ccsBase + amendCount++;
					}
					featIdx0 = featIdx1;
				}
			}
			rDomAligns.push_back(alignRec);
		}
	}
	
	// -- added 9/9/2014 handling structure motifs -- attach to 
	for (size_t iidx = 0, iidxEnd = m_vecSDIndice.size(); iidx < iidxEnd; ++iidx)
	{
		__TOflAlignRecord rec;
		
		rec.idx = m_vecSDIndice[iidx];
		
		rec.pAlign =  &(rAlignments[rec.idx]);
		
		rec.pCDInfo = rDomInfo.FindCDInfo(rec.pAlign->m_uiPSSMID);
		
		rec.pClst = NULL;
		
		rec.pRootCDInfo = NULL;
		
		rec.idxidx = ccsBase + amendCount++;
		
		rFeatAligns.push_back(rec);
	}
}


void TOflAlignIndice::ExtractFeatAligns(const __TOflAlignRecord &rRepRec, const vector<TOflCDAlignInfo> &rAlignments, const COflDomClstInfo & rDomInfo, vector<__TOflAlignRecord> &rResult) const
{


	__TOflAlignRecord rec;
	if (rRepRec.pCDInfo->m_uiClusterPSSMID > 0)	//has cluster, do cluster match
	{
		for (size_t fiidx = 0, fiidxEnd = m_vecQualifiedFeatIndice.size(); fiidx < fiidxEnd; ++fiidx)
		{
			rec.pAlign = &(rAlignments[m_vecQualifiedFeatIndice[fiidx]]);
			if (rec.pAlign->m_iRegionIdx == rRepRec.pAlign->m_iRegionIdx)	//region match
			{
				rec.pCDInfo = rDomInfo.FindCDInfo(rec.pAlign->m_uiPSSMID);
				if (rec.pCDInfo->m_uiClusterPSSMID == rRepRec.pCDInfo->m_uiClusterPSSMID)	//matched cluster
				{
					rec.pClst = rRepRec.pClst;
					rec.idx = m_vecQualifiedFeatIndice[fiidx];
					//rec.idxidx = FeatIdx2Iidx(rec.idx);
					if (rec.pCDInfo->m_uiHierarchyRoot > 0)
					{
						if (rec.pCDInfo->m_uiHierarchyRoot == rec.pAlign->m_uiPSSMID) rec.pRootCDInfo = rec.pCDInfo;
						else rec.pRootCDInfo = rDomInfo.FindCDInfo(rec.pCDInfo->m_uiHierarchyRoot);
					}
					rResult.push_back(rec);
				}
			}
		}
	}
	else	//no cluster, just match region class -- should not happen
	{
		for (size_t fiidx = 0, fiidxEnd = m_vecQualifiedFeatIndice.size(); fiidx < fiidxEnd; ++fiidx)
		{
			rec.pAlign = &(rAlignments[m_vecQualifiedFeatIndice[fiidx]]);
			if (rec.pAlign->m_iRegionIdx == rRepRec.pAlign->m_iRegionIdx)	//region match
			{
				rec.pCDInfo = rDomInfo.FindCDInfo(rec.pAlign->m_uiPSSMID);
				rec.pClst = rDomInfo.FindClusterInfo(rec.pCDInfo->m_uiClusterPSSMID);
				if (rec.pCDInfo->m_uiHierarchyRoot > 0)
				{
					if (rec.pCDInfo->m_uiHierarchyRoot == rec.pAlign->m_uiPSSMID) rec.pRootCDInfo = rec.pCDInfo;
					else rec.pRootCDInfo = rDomInfo.FindCDInfo(rec.pCDInfo->m_uiHierarchyRoot);
				}
				rec.idx = m_vecQualifiedFeatIndice[fiidx];
				//rec.idxidx = FeatIdx2Iidx(rec.idx);
				rResult.push_back(rec);
			}
		}
	}
}

/**********************************************************************
*	Biodata structure -- Sequence information
***********************************************************************/
struct TOflSeqInfo
{
	int m_iGi;
	string m_strQueryID;
	unsigned int m_uiSeqLen;	//length of original sequence, could be na
	bool m_bIsProtein;
	string m_strDefline;
	string m_strMessage;
	
	string m_dimTranslated[TOTAL_RFS];
	
	TOflSeqInfo(void): m_iGi(0), m_strQueryID(k_strEmptyString), m_uiSeqLen(0), m_bIsProtein(true), m_strDefline(k_strEmptyString), m_strMessage(k_strEmptyString) {};

	void Get1LtrSeq(string &dst) const;
	size_t Translate(const TOflAlignInfo& alignInfo, string &dst) const;
	unsigned int GetSeqLength(void) const {return m_uiSeqLen;}
	
	// -- internals
	void NACommit(void);
	void PRCommit(void);
};

void TOflSeqInfo::Get1LtrSeq(string &dst) const
{
	dst.clear();
	if (m_bIsProtein) dst = m_dimTranslated[0];
	else dst.append(m_uiSeqLen, 'N');
}

size_t TOflSeqInfo::Translate(const TOflAlignInfo& alignInfo, string &dst) const
{
	dst = m_dimTranslated[alignInfo.GetRFIdx()];
	return m_uiSeqLen;
}

void TOflSeqInfo::NACommit(void)
{
	m_bIsProtein = false;
	size_t aaLen = m_uiSeqLen / 3;
	
	m_dimTranslated[0].clear();
	m_dimTranslated[0].append(aaLen, '-');
	
	m_dimTranslated[3] = m_dimTranslated[0];
	
	aaLen = (m_uiSeqLen - 1) / 3;
	
	m_dimTranslated[1].clear();
	m_dimTranslated[1].append(aaLen, '-');
	
	m_dimTranslated[4] = m_dimTranslated[1];
	
	aaLen = (m_uiSeqLen - 2) / 3;
	
	m_dimTranslated[2].clear();
	m_dimTranslated[2].append(aaLen, '-');
	m_dimTranslated[5] = m_dimTranslated[2];
}

void TOflSeqInfo::PRCommit(void)
{
	m_bIsProtein = true;
	m_dimTranslated[0].clear();
	m_dimTranslated[0].append(m_uiSeqLen, '-');
}

/**********************************************************************
*	Biodata structure -- CD Query result - protein and na
***********************************************************************/
struct TOflCDQuery: public TOflSeqInfo, public TOflAlignIndice
{
	vector<TOflCDAlignInfo> m_vecAlignments;
	
	TOflCDQuery(void);
	
	void CreateRecordSets(const COflDomClstInfo & rDomInfo, vector<__TOflAlignRecord> &rDomAligns, vector<__TOflAlignRecord> &rFeatAligns, int mode) const;
	
	// -- calculate segments from 
	void ParseAlignStrings(const string &qseq, const string &hseq, int qf, int hf, size_t start, size_t n, TOflCDAlignInfo &dst, int rfidx = 0);
	
};



// -- extend for nucleic acid queries
struct TOflCdQueryEx: public TOflCDQuery
{
	TOflAlignIndice m_dimSplitAligns[TOTAL_RFS];
	TOflCdQueryEx(void): TOflCDQuery(), m_dimSplitAligns() {};
	
	//void CreateRecordSetsEx(int rfidx, const COflDomClstInfo & rDomInfo, vector<__TOflAlignRecord> &rDomAligns, vector<__TOflAlignRecord> &rFeatAligns, bool bConcise = true) const;
	void CreateRecordSetsEx(int rfidx, const COflDomClstInfo & rDomInfo, vector<__TOflAlignRecord> &rDomAligns, vector<__TOflAlignRecord> &rFeatAligns, int mode) const;
	size_t CollectCdQueries(list<TOflCDQuery> &dst) const;
		
	void ProcessQueryResult(CRef< CIteration > refOflQueryResult, double evcut = 0.01);
		
	void Print(int idxBlObj, const COflDomClstInfo &domInfo, ostream &os, int dmode = TDataModes::e_rep) const;
};

TOflCDQuery::TOflCDQuery(void):
	TOflSeqInfo(), TOflAlignIndice(), m_vecAlignments()
{}

void TOflCDQuery::CreateRecordSets(const COflDomClstInfo & rDomInfo, vector<TOflAlignIndice::__TOflAlignRecord> &rDomAligns, vector<TOflAlignIndice::__TOflAlignRecord> &rFeatAligns, int mode) const
{
	TOflAlignIndice::CreateRecordSets(m_vecAlignments, rDomInfo, rDomAligns, rFeatAligns, mode);
}

void TOflCDQuery::ParseAlignStrings(const string &qseq, const string &hseq, int qf, int hf, size_t start, size_t n, TOflCDAlignInfo &dst, int rfidx)
{
	unsigned int segLen = 0;
	bool bEffSeg = true;

	while (start < n)
	{

		if (bEffSeg)	//initially, both chains are in effect
		{
			if ('-' != qseq[start] && '-' != hseq[start])
				++segLen;
			else
			{

				bEffSeg = false;
				
				dst.m_vecMStarts.push_back(qf);
				dst.m_vecSStarts.push_back(hf);
				dst.m_vecLens.push_back(segLen);
				
				// -- write in master sequence data
				m_dimTranslated[rfidx].replace(qf, segLen, qseq, start - segLen, segLen);
				
				qf += segLen;
				hf += segLen;
				segLen = 0;
				if ('-' == qseq[start])
					--qf;
				else
				{
					--hf;
					//write in one more residue
					m_dimTranslated[rfidx][qf] = qseq[start];
				}
			}
		}
		else	//not in effective seg
		{
			//!bEffSeg
			bool segStart = true;
			if ('-' != qseq[start])
			{
				++qf;
				m_dimTranslated[rfidx][qf] = qseq[start];
			}
			else
				segStart = false;
				
			if ('-' != hseq[start])
				++hf;
			else
				segStart = false;
				
			if ((bEffSeg = segStart)) ++segLen;
			
		}
			
		++start;
	}//while
	// -- last segment
	if (bEffSeg && segLen > 0)
	{

		dst.m_vecMStarts.push_back(qf);
		dst.m_vecSStarts.push_back(hf);
		dst.m_vecLens.push_back(segLen);
		
		m_dimTranslated[rfidx].replace(qf, segLen, qseq, start - segLen, segLen);
	}
}

void TOflCdQueryEx::CreateRecordSetsEx(int rfidx, const COflDomClstInfo & rDomInfo, vector<__TOflAlignRecord> &rDomAligns, vector<__TOflAlignRecord> &rFeatAligns, int mode) const
{
	m_dimSplitAligns[rfidx].CreateRecordSets(m_vecAlignments, rDomInfo, rDomAligns, rFeatAligns, mode);
}

size_t TOflCdQueryEx::CollectCdQueries(list<TOflCDQuery> &dst) const
{
	dst.clear();
	size_t total = 0;
	if (m_bIsProtein)
	{
		dst.push_back((const TOflCDQuery &)(*this));
		return 1;
	}
	else	//push as a serious of TOflCDQuery
	{
		for (int rf = 0; rf < TOTAL_RFS; ++rf)
		{
			if (!m_dimSplitAligns[rf].m_vecSortedIndice.empty())
			{
				list<TOflCDQuery> :: iterator iterNewRf = dst.insert(dst.end(), (const TOflCDQuery &)(*this));
				++total;
				iterNewRf->m_vecSortedIndice = m_dimSplitAligns[rf].m_vecSortedIndice;
				iterNewRf->m_vecConciseIndice = m_dimSplitAligns[rf].m_vecConciseIndice;
				iterNewRf->m_vecStdIndice = m_dimSplitAligns[rf].m_vecStdIndice;
				iterNewRf->m_strDefline = RF_TITLES[rf];
				// -- add reading frame id
			}
		}
	}
	return total;
}


void TOflCdQueryEx::ProcessQueryResult(CRef<CIteration> refOflQueryResult, double evcut)
{
	m_vecAlignments.clear();
	
	if (refOflQueryResult->CanGetQuery_ID())
		m_strQueryID = refOflQueryResult->GetQuery_ID();
	
	if (refOflQueryResult->CanGetQuery_def())
		m_strDefline = refOflQueryResult->GetQuery_def();	//length of original sequence, could be na
	
	if (!g_bSilent) cerr << "Processing Query " << m_strDefline << endl;
	
	if (refOflQueryResult->CanGetQuery_len())
		m_uiSeqLen = refOflQueryResult->GetQuery_len();	//length of original sequence, could be na
		
	
	if (refOflQueryResult->CanGetMessage())
		m_strMessage = refOflQueryResult->GetMessage();	//length of original sequence, could be na
		
	if (refOflQueryResult->CanGetHits())
	{
		const CIteration::THits &hits = refOflQueryResult->GetHits();	//list< CRef< CHit > >
		if (!hits.empty())
		{


			m_vecAlignments.reserve(hits.size() * 3);	//estimated
			
			TOflCDAlignInfo alnVal;
			
			for (CIteration::THits::const_iterator iterHit = hits.begin(), iterHitEnd = hits.end(); iterHitEnd != iterHit; ++iterHit)
			{
				unsigned int pssmid = (unsigned int)atoi((*iterHit)->GetAccession().c_str());
				
				unsigned int cdLen = (*iterHit)->GetLen();
				const CHit::THsps &hsps = (*iterHit)->GetHsps();	//list< CRef< CHsp > > THsps
				if (!g_bSilent) cerr << "Domain hit " << pssmid << ", total occurrance = " << hsps.size() << endl;	
				for (CHit::THsps::const_iterator iterHsp = hsps.begin(), iterHspEnd = hsps.end(); iterHspEnd != iterHsp; ++iterHsp)
				{
					double eval = (*iterHsp)->GetEvalue();
					if (eval > evcut) continue;
					
					m_vecAlignments.push_back(alnVal);
					TOflCDAlignInfo &dst = *m_vecAlignments.rbegin();
					
					dst.m_uiPSSMID = pssmid;
					dst.m_uiAlignedLen = (*iterHsp)->GetAlign_len();
					dst.m_iScore = (int)((*iterHsp)->GetScore() + 0.5);
					dst.m_dEValue = eval;
					dst.m_dBitScore = (*iterHsp)->GetBit_score();
					dst.m_iNumIdent = (*iterHsp)->GetIdentity();
					dst.m_dSeqIdentity = (double)dst.m_iNumIdent / (double)cdLen * 100.0;
					dst.m_iReadingFrame = (*iterHsp)->GetQuery_frame();
					
					if (dst.m_iReadingFrame > 1 || dst.m_iReadingFrame < 0)	//na determined
					{
						if (m_dimTranslated[0].empty())	//not filled yet
							NACommit();
					}
					
					
					dst.m_iFrom = (*iterHsp)->GetQuery_from() - COORDSBASE;
					dst.m_iTo = (*iterHsp)->GetQuery_to() - COORDSBASE;

					// -- start to parse sequence
					const string &qseq = (*iterHsp)->GetQseq(), hseq = (*iterHsp)->GetHseq();
					size_t alnLen = qseq.size(), alnResIdx = 0;	//assume qseq.size() == hseq.size()

					// -- incase of blast problem...
					while (alnResIdx < alnLen && ('-' == qseq[alnResIdx] || '-' == hseq[alnResIdx])) ++alnResIdx;
					
					// -- see if NA
					if (m_dimTranslated[0].empty())	//not filled yet
					{
						int aaRes = 0, t0 = alnResIdx;
						
						while (t0 < alnLen)
						{
							if ('-' == qseq[t0])
								++t0;
							else
							{
								size_t t1 = qseq.find('-', t0);
								if (string::npos == t1)	//
									t1 = alnLen;
								aaRes += t1 - t0;
								t0 = t1 + 1;
							}
						}
						if (aaRes == dst.m_iTo - dst.m_iFrom + 1)
							PRCommit();
						else if (aaRes * 3 == dst.m_iTo - dst.m_iFrom + 1)
							NACommit();
					}
					if (m_bIsProtein)
					{
						dst.m_iReadingFrame = 0;
						dst.m_eAlignType = TOflAlignInfo::eNormal;
						dst.m_eStrand = eNa_strand_unknown;
						int qfrom = dst.m_iFrom, hfrom = (*iterHsp)->GetHit_from() - COORDSBASE;
						
						ParseAlignStrings(qseq, hseq, qfrom, hfrom, alnResIdx, alnLen, dst, 0);
					}
					else	//na
					{
						// -- calculate 
						int rfidx = RF2Idx(dst.m_iReadingFrame);
						dst.m_eAlignType = TOflAlignInfo::ePr2Na;

						if (dst.m_iReadingFrame > 0)	//plus
						{
							dst.m_eStrand = eNa_strand_plus;
							dst.m_iReadingFrame -= 1;	//from 1, 2, 3 to 0, 1 ,2
							int qfrom = dst.NAPlus2Pr(dst.m_iFrom), hfrom = (*iterHsp)->GetHit_from() - COORDSBASE;

							ParseAlignStrings(qseq, hseq, qfrom, hfrom, alnResIdx, alnLen, dst, rfidx);

						}	//plus strand
						else	//minus strand
						{
							dst.m_eStrand = eNa_strand_minus;
							// -- squeeze m_uiSeqLen in dst.m_iReadingFrame, rf -1, -2, -3 -> 0, 1, 2
							dst.m_iReadingFrame = (-dst.m_iReadingFrame - 1) | (m_uiSeqLen << 2);
							int qfrom = dst.NAMinus2Pr(dst.m_iTo), hfrom = (*iterHsp)->GetHit_from() - COORDSBASE;

							ParseAlignStrings(qseq, hseq, qfrom, hfrom, alnResIdx, alnLen, dst, rfidx);
						}
					}
				}
			}
		}
	}
}

// -- helper function
void PrintDomLine(ostream &os, const char *pType, const TOflAlignIndice::__TOflAlignRecord &rec)
{
	os << pType << DELIMIT << rec.pAlign->m_uiPSSMID << DELIMIT << rec.pAlign->m_iFrom + COORDSBASE << DELIMIT << rec.pAlign->m_iTo + COORDSBASE << DELIMIT << rec.pAlign->m_dEValue << DELIMIT << rec.pAlign->m_dBitScore << DELIMIT << rec.pCDInfo->m_strAccession << DELIMIT << rec.pCDInfo->m_strShortName << DELIMIT;
	bool bTmMiss = false;
	if (rec.pAlign->m_dNMissing > 0.2)
	{
		os << 'N';
		bTmMiss = true;
	}
	if (rec.pAlign->m_dCMissing > 0.2)
	{
		os << 'C';
		bTmMiss = true;
	}
	if (!bTmMiss)
		os << '-';
	os << DELIMIT << (NULL == rec.pClst ? '-' :  rec.pCDInfo->m_uiClusterPSSMID);
}

void PrintClstLine(ostream &os, const TOflAlignIndice::__TOflAlignRecord &rec)
{
	os << HITTYPE_CLUSTER << DELIMIT << rec.pCDInfo->m_uiClusterPSSMID << DELIMIT << rec.pAlign->m_iFrom + COORDSBASE << DELIMIT << rec.pAlign->m_iTo + COORDSBASE << DELIMIT << rec.pAlign->m_dEValue << DELIMIT << rec.pAlign->m_dBitScore << DELIMIT << rec.pClst->m_strAccession << DELIMIT << rec.pClst->m_strShortName << DELIMIT;
	bool bTmMiss = false;
	if (rec.pAlign->m_dNMissing > 0.2)
	{
		os << 'N';
		bTmMiss = true;
	}
	if (rec.pAlign->m_dCMissing > 0.2)
	{
		os << 'C';
		bTmMiss = true;
	}
	if (!bTmMiss)
		os << '-';
	os << DELIMIT << '-';
}

struct __MotifType
{
	list<TOflCDFeat> :: const_iterator iterMotifFeat;
	vector<TOflAlignIndice::__TOflAlignRecord> :: const_iterator iterAlignRec;
	bool bIsSpecific;
	unsigned int uiSrcPSSMId;
	__MotifType(list<TOflCDFeat> :: const_iterator iterM, vector<TOflAlignIndice::__TOflAlignRecord> :: const_iterator iterA, bool spec, unsigned int srcpssm): 
		iterMotifFeat(iterM), iterAlignRec(iterA), bIsSpecific(spec), uiSrcPSSMId(srcpssm) {};
};


void TOflCdQueryEx::Print(int idxBlObj, const COflDomClstInfo &domInfo, ostream &os, int dmode) const
{
	os << QUERYSTART << DELIMIT << m_strQueryID << DELIMIT;
//const char * const QUERY_TYPE_PEPTIDE = "Peptide";
//const char * const QUERY_TYPE_NUCLEOTIDE = "Nucleotide";
	if (m_bIsProtein)
	{
		os << QUERY_TYPE_PEPTIDE << DELIMIT << m_uiSeqLen << DELIMIT << m_strDefline << endl;
		vector<TOflAlignIndice::__TOflAlignRecord> vecDomAligns, vecFeatAligns;	//vecDomFams for output of superfamilies of domain hits
		CreateRecordSets(domInfo, vecDomAligns, vecFeatAligns, dmode);
		
		map<int, TOflAlignIndice::__TOflAlignRecord > mapDomFams;


		if (TTargetData::e_feats != g_iTargetData)
		{
			if (!vecDomAligns.empty())
			{
				os << DOMSTART << endl;
				if (TDataModes::e_rep == dmode)
				{
					for (vector<TOflAlignIndice::__TOflAlignRecord> :: const_iterator iterRec = vecDomAligns.begin(), iterRecEnd = vecDomAligns.end(); iterRecEnd != iterRec; ++iterRec)
					{
						if (NULL == iterRec->pCDInfo) continue;
						os << idxBlObj << DELIMIT << m_strQueryID << DELIMIT;
						if (1 == iterRec->pAlign->m_iRepClass)	//multi, should be no more
							PrintDomLine(os, HITTYPE_MULTIDOM, *iterRec);
						else if (iterRec->pAlign->m_bSpecQualified)
						{
							PrintDomLine(os, HITTYPE_SPECIFIC, *iterRec);
							if (NULL != iterRec->pClst)
							{
								int key = iterRec->pAlign->m_iRegionIdx * OVERLAPLEADING + iterRec->pCDInfo->m_uiClusterPSSMID;
								if (mapDomFams.end() == mapDomFams.find(key))
									mapDomFams.emplace(key, *iterRec);
							}
						}
						else if (NULL != iterRec->pClst)
							PrintClstLine(os, *iterRec);
						else	//not qualified for specific, but no cluster -- non-specific
							PrintDomLine(os, HITTYPE_NONSPECIFIC, *iterRec);
						os << endl;
					}	//dom loops end
				}	//rep mode
				else	//nonrep mode
				{
					for (vector<TOflAlignIndice::__TOflAlignRecord> :: const_iterator iterRec = vecDomAligns.begin(), iterRecEnd = vecDomAligns.end(); iterRecEnd != iterRec; ++iterRec)
					{
						if (NULL == iterRec->pCDInfo) continue;
						os << idxBlObj << DELIMIT << m_strQueryID << DELIMIT;
						if (1 == iterRec->pAlign->m_iRepClass)	//multi
							PrintDomLine(os, HITTYPE_MULTIDOM, *iterRec);
						else if (iterRec->pAlign->m_bRep && iterRec->pAlign->m_bSpecQualified)
							PrintDomLine(os, HITTYPE_SPECIFIC, *iterRec);
						else
							PrintDomLine(os, HITTYPE_NONSPECIFIC, *iterRec);
						if (NULL != iterRec->pClst)
						{
							int key = iterRec->pAlign->m_iRegionIdx * OVERLAPLEADING + iterRec->pCDInfo->m_uiClusterPSSMID;
							if (mapDomFams.end() == mapDomFams.find(key))
								mapDomFams.emplace(key, *iterRec);
						}
						os << endl;
					}
				}
				
				os << DOMEND << endl;
				
				if (g_bSuperfams)	//user want to see all superfams
				{
					os << FAMSTART << endl;
					for (const auto & v : mapDomFams)
					{
						os << idxBlObj << DELIMIT << m_strQueryID << DELIMIT;
						PrintClstLine(os, v.second);
						os << endl;
					}
					os << FAMEND << endl;
				}
			}
		}
		
		if (TTargetData::e_doms != g_iTargetData)
		{
			vector<__MotifType> vecMotifRecs;
			
			vector<TOflCDQuery::__TOflAlignRecord> :: const_iterator iterFeatAlign = vecFeatAligns.begin(), iterFeatAlignEnd = vecFeatAligns.end();
			
			if (iterFeatAlignEnd != iterFeatAlign)	//not empty())
			{
				if (!iterFeatAlign->pCDInfo->m_bIsStructDom)	//has regular domains
				{
					os << FEATSTART << endl;
					for ( ; iterFeatAlignEnd != iterFeatAlign; ++iterFeatAlign)
					{
						if (iterFeatAlign->pCDInfo->m_bIsStructDom) break;
							
						bool bIsSpecific = iterFeatAlign->pAlign->m_bRep && iterFeatAlign->pAlign->m_bSpecQualified;
						const char * pType = ANNOTTYPE_SPECIFIC;
						const list <TOflCDFeat> * pFeats = &(iterFeatAlign->pCDInfo->m_lstSpecFeatures);;
						unsigned int uiSrcPSSMId = iterFeatAlign->pAlign->m_uiPSSMID;
						
						if (!bIsSpecific && NULL != iterFeatAlign->pRootCDInfo)
						{
							pType = ANNOTTYPE_GENERIC;
							pFeats = &(iterFeatAlign->pCDInfo->m_lstGenFeatures);
							uiSrcPSSMId = iterFeatAlign->pCDInfo->m_uiHierarchyRoot;
						}
							
						for (list<TOflCDFeat> :: const_iterator iterFeat = pFeats->begin(); iterFeat != pFeats->end(); ++iterFeat)
						{
							if (TOflCDFeat::eType_StructMotif == iterFeat->m_iType)
							{
								vecMotifRecs.push_back(__MotifType(iterFeat, iterFeatAlign, bIsSpecific, uiSrcPSSMId));
								continue;
							}
							
							CSegSet featsegs(*iterFeat);
							iterFeatAlign->pAlign->MapSegSet(featsegs, false);
							
							if (!featsegs.IsEmpty() && ((double)featsegs.GetTotalResidues() / (double)iterFeat->GetCompleteSize() > 0.8))
							{
								vector<CSegSet::TResiduePos> vecMappedPos;
								int rfidx = iterFeatAlign->pAlign->GetTranslatedPosMap(featsegs, vecMappedPos);
    	        	
								if (iterFeat->MotifCheck(vecMappedPos, m_dimTranslated[rfidx]) > 0) continue;	//failed motif check
    	        	
									
								os << idxBlObj << DELIMIT << m_strQueryID << DELIMIT << pType << DELIMIT << iterFeat->m_strTitle << DELIMIT;
								
								const TSegs& segs = featsegs.GetSegs();
								
								size_t mapped = 0;
								
								char dimDelimit[2] = {0, 0};
								for (TSegs::const_iterator iterSeg = segs.begin(); iterSeg != segs.end(); ++iterSeg)
								{
									for (int res = iterSeg->from; res <= iterSeg->to; ++res)
									{
										os << dimDelimit << (char)toupper(m_dimTranslated[rfidx][res]) << (res + COORDSBASE);
										dimDelimit[0] = COORDELIMIT;
										++mapped;
									}
								}
								os << DELIMIT << iterFeat->GetTotalResidues() << DELIMIT << mapped << DELIMIT << uiSrcPSSMId << endl;
							}
						}
					}
					
					// -- non-motif features from SD
					
					for ( ; iterFeatAlignEnd != iterFeatAlign; ++iterFeatAlign)
					{
						//if (iterFeatAlign->pCDInfo->m_bIsStructDom) break;
							
						//bool bIsSpecific = iterFeatAlign->pAlign->m_bRep && iterFeatAlign->pAlign->m_bSpecQualified;
						//const char * pType = ANNOTTYPE_SPECIFIC;
						//const list <TOflCDFeat> * pFeats = &(iterFeatAlign->pCDInfo->m_lstSpecFeatures);;
						//unsigned int uiSrcPSSMId = iterFeatAlign->pAlign->m_uiPSSMID;
						
						//if (!bIsSpecific && NULL != iterFeatAlign->pRootCDInfo)
						//{
						//	pType = ANNOTTYPE_GENERIC;
						//	pFeats = &(iterFeatAlign->pCDInfo->m_lstGenFeatures);
						//	uiSrcPSSMId = iterFeatAlign->pCDInfo->m_uiHierarchyRoot;
						//}
							
						for (list<TOflCDFeat> :: const_iterator iterFeat = iterFeatAlign->pCDInfo->m_lstSpecFeatures.begin(); iterFeat != iterFeatAlign->pCDInfo->m_lstSpecFeatures.end(); ++iterFeat)
						{
							if (TOflCDFeat::eType_StructMotif == iterFeat->m_iType)
							{
								vecMotifRecs.push_back(__MotifType(iterFeat, iterFeatAlign, true, 0));
								continue;
							}
							
							CSegSet featsegs(*iterFeat);
							iterFeatAlign->pAlign->MapSegSet(featsegs, false);
							
							if (!featsegs.IsEmpty() && ((double)featsegs.GetTotalResidues() / (double)iterFeat->GetCompleteSize() > 0.8))
							{
								vector<CSegSet::TResiduePos> vecMappedPos;
								int rfidx = iterFeatAlign->pAlign->GetTranslatedPosMap(featsegs, vecMappedPos);
    	        	
								if (iterFeat->MotifCheck(vecMappedPos, m_dimTranslated[rfidx]) > 0) continue;	//failed motif check
    	        	
									
								os << idxBlObj << DELIMIT << m_strQueryID << DELIMIT << ANNOTTYPE_SPECIFIC << DELIMIT << iterFeat->m_strTitle << DELIMIT;
								
								const TSegs& segs = featsegs.GetSegs();
								
								size_t mapped = 0;
								
								char dimDelimit[2] = {0, 0};
								for (TSegs::const_iterator iterSeg = segs.begin(); iterSeg != segs.end(); ++iterSeg)
								{
									for (int res = iterSeg->from; res <= iterSeg->to; ++res)
									{
										os << dimDelimit << (char)toupper(m_dimTranslated[rfidx][res]) << (res + COORDSBASE);
										dimDelimit[0] = COORDELIMIT;
										++mapped;
									}
								}
								os << DELIMIT << iterFeat->GetTotalResidues() << DELIMIT << mapped << DELIMIT << 0/*iterFeatAlign->pAlign->m_uiPSSMID*/ << endl;
							}
						}
					}
					
					os << FEATEND << endl;
				}
				
				// -- Motif part
				if (!vecMotifRecs.empty())
				{
					// -- start SD
					// -- features
					os << MOTIFSTART << endl;
					vector<__MotifType> :: const_iterator iterM = vecMotifRecs.begin(), iterMEnd = vecMotifRecs.end();
					for ( ; iterMEnd != iterM; ++iterM)
					{
						if (iterM->iterAlignRec->pCDInfo->m_bIsStructDom) break;
							
						CSegSet featsegs(*iterM->iterMotifFeat);
						iterM->iterAlignRec->pAlign->MapSegSet(featsegs);
						
						if (!featsegs.IsEmpty() && ((double)featsegs.GetTotalResidues() / (double)iterM->iterMotifFeat->GetCompleteSize() > 0.8))
						{
							os << idxBlObj << DELIMIT << m_strQueryID << DELIMIT << iterM->iterMotifFeat->m_strTitle << DELIMIT;
							int res0 = featsegs.GetLeft(), res1 = featsegs.GetRight(), rfidx = iterM->iterAlignRec->pAlign->GetRFIdx();
							os << (char)toupper(m_dimTranslated[rfidx][res0]) << res0 + COORDSBASE << DELIMIT << res1 + COORDSBASE << (char)toupper(m_dimTranslated[rfidx][res1]) << DELIMIT << iterM->uiSrcPSSMId << endl;
						}
					}
					
					// Motif from SD
					for ( ; iterMEnd != iterM; ++iterM)
					{
						CSegSet featsegs(*iterM->iterMotifFeat);
						if (!iterM->iterAlignRec->pAlign->m_ClipSet.IsEmpty())
							featsegs.Cross(iterM->iterAlignRec->pAlign->m_ClipSet);
						
						iterM->iterAlignRec->pAlign->MapSegSet(featsegs);
						
						if (!featsegs.IsEmpty() && ((double)featsegs.GetTotalResidues() / (double)iterM->iterMotifFeat->GetCompleteSize() > 0.8))
						{
							os << idxBlObj << DELIMIT << m_strQueryID << DELIMIT << iterM->iterMotifFeat->m_strTitle << DELIMIT;
							int res0 = featsegs.GetLeft(), res1 = featsegs.GetRight(), rfidx = iterM->iterAlignRec->pAlign->GetRFIdx();
							os << (char)toupper(m_dimTranslated[rfidx][res0]) << res0 + COORDSBASE << DELIMIT << res1 + COORDSBASE << (char)toupper(m_dimTranslated[rfidx][res1]) << DELIMIT << iterM->uiSrcPSSMId << endl;
						}
					}
					
					// -- now start SD
					//for ( ; iterFeatAlignEnd != iterFeatAlign; ++iterFeatAlign)
					//{
					//	int rfidx = iterFeatAlign->pAlign->GetRFIdx();
					//	
					//	for (list<TOflCDFeat> :: const_iterator iterFeat = iterFeatAlign->pCDInfo->m_lstSpecFeatures.begin(); iterFeat != iterFeatAlign->pCDInfo->m_lstSpecFeatures.end(); ++iterFeat)
					//	{
					//		
					//		CSegSet featsegs(*iterFeat);
					//		iterFeatAlign->pAlign->MapSegSet(featsegs);
					//		if (!featsegs.IsEmpty())
					//		{
					//			os << idxBlObj << DELIMIT << m_strQueryID << DELIMIT << iterFeat->m_strTitle << DELIMIT;
					//			//int rfidx = iterFeatAlign->pAlign->GetRFIdx();
					//			int res0 = featsegs.GetLeft(), res1 = featsegs.GetRight();
					//			os << (char)toupper(m_dimTranslated[rfidx][res0]) << res0 + COORDSBASE << DELIMIT << res1 + COORDSBASE << (char)toupper(m_dimTranslated[rfidx][res1]) << DELIMIT << '-' << endl;
					//		}
					//	}
					//}
					
					os << MOTIFEND << endl;
				}
			}
		}
	}
	else	//na
	{
		os << QUERY_TYPE_NUCLEOTIDE << DELIMIT << m_uiSeqLen << DELIMIT << m_strDefline << endl;
		
		vector<TOflAlignIndice::__TOflAlignRecord> dimDomAligns[TOTAL_RFS], dimFeatVecs[TOTAL_RFS];
		
		bool bHasFeats = false, bHasRegFeats = false;

		for (int rfidx = 0; rfidx < TOTAL_RFS; ++rfidx)
		{
			CreateRecordSetsEx(rfidx, domInfo, dimDomAligns[rfidx], dimFeatVecs[rfidx], dmode);
			bHasFeats = bHasFeats || !dimFeatVecs[rfidx].empty();
			bHasRegFeats = (bHasRegFeats || !(dimFeatVecs[rfidx].empty() || dimFeatVecs[rfidx][0].pCDInfo->m_bIsStructDom));
		}

		if (TTargetData::e_feats != g_iTargetData)
		{
			map<int, TOflAlignIndice::__TOflAlignRecord > dimDomFams[TOTAL_RFS];
			
			os << DOMSTART << endl;
			for (int rfidx = 0; rfidx < TOTAL_RFS; ++rfidx)
			{
				int rf = Idx2RF(rfidx);
				if (!dimDomAligns[rfidx].empty())
				{
					if (TDataModes::e_rep == dmode)
					{
						for (vector<TOflAlignIndice::__TOflAlignRecord> :: const_iterator iterRec = dimDomAligns[rfidx].begin(), iterRecEnd = dimDomAligns[rfidx].end(); iterRecEnd != iterRec; ++iterRec)
						{
							if (NULL == iterRec->pCDInfo) continue;
							os << idxBlObj << DELIMIT << m_strQueryID << '[' << rf << ']' << DELIMIT;
							if (1 == iterRec->pAlign->m_iRepClass)	//multi
								PrintDomLine(os, HITTYPE_MULTIDOM, *iterRec);
							else if (iterRec->pAlign->m_bSpecQualified)
							{
								PrintDomLine(os, HITTYPE_SPECIFIC, *iterRec);
								if (NULL != iterRec->pClst)
								{
									int key = iterRec->pAlign->m_iRegionIdx * OVERLAPLEADING + iterRec->pCDInfo->m_uiClusterPSSMID;
									if (dimDomFams[rfidx].end() == dimDomFams[rfidx].find(key))
										dimDomFams[rfidx].emplace(key, *iterRec);
								}
							}
							else if (NULL != iterRec->pClst)
								PrintClstLine(os, *iterRec);
							else	//not qualified for specific, but no cluster -- non-specific
								PrintDomLine(os, HITTYPE_NONSPECIFIC, *iterRec);
							os << endl;
						}	//dom loops end
					}	//rep mode
					else	//nonrep mode
					{
						for (vector<TOflAlignIndice::__TOflAlignRecord> :: const_iterator iterRec = dimDomAligns[rfidx].begin(), iterRecEnd = dimDomAligns[rfidx].end(); iterRecEnd != iterRec; ++iterRec)
						{
							if (NULL == iterRec->pCDInfo) continue;
							os << idxBlObj << DELIMIT << m_strQueryID << '[' << rf << ']' << DELIMIT;
							if (1 == iterRec->pAlign->m_iRepClass)	//multi
								PrintDomLine(os, HITTYPE_MULTIDOM, *iterRec);
							else if (iterRec->pAlign->m_bRep && iterRec->pAlign->m_bSpecQualified)
								PrintDomLine(os, HITTYPE_SPECIFIC, *iterRec);
							else
								PrintDomLine(os, HITTYPE_NONSPECIFIC, *iterRec);
								
							if (NULL != iterRec->pClst)
							{
								int key = iterRec->pAlign->m_iRegionIdx * OVERLAPLEADING + iterRec->pCDInfo->m_uiClusterPSSMID;
								if (dimDomFams[rfidx].end() == dimDomFams[rfidx].find(key))
									dimDomFams[rfidx].emplace(key, *iterRec);
							}
							
							os << endl;
						}
					}
				}
			}
			os << DOMEND << endl;
			
			if (g_bSuperfams)	//user want to see all superfams
			{
				os << FAMSTART << endl;
				for (int rfidx = 0; rfidx < TOTAL_RFS; ++rfidx)
				{
					int rf = Idx2RF(rfidx);
					
					for (const auto & v : dimDomFams[rfidx])
					{
						os << idxBlObj << DELIMIT << m_strQueryID << DELIMIT;
						PrintClstLine(os, v.second);
						os << endl;
					}
					
				}
				os << FAMEND << endl;
			}
		}
			
		if (TTargetData::e_doms != g_iTargetData)
		{
			if (bHasFeats)
			{
				vector<__MotifType> dimMotifTypes[TOTAL_RFS];
//				vector<TOflCDQuery::__TOflAlignRecord> :: const_iterator dimSDIter[TOTAL_RFS];
				
				bool bHasMotifs = false;
				
				if (bHasRegFeats)
					os << FEATSTART << endl;
				
				int iNegRF = TOTAL_RFS / 2;
				for (int rfidx = 0; rfidx < iNegRF; ++rfidx)	//plus strand reading frames
				{
					int rf = rfidx + 1;
					vector<TOflCDQuery::__TOflAlignRecord> :: const_iterator iterFeatAlign = dimFeatVecs[rfidx].begin(), iterFeatAlignEnd = dimFeatVecs[rfidx].end();
					//dimSDIter[rfidx] = iterFeatAlignEnd;
					for ( ; iterFeatAlignEnd != iterFeatAlign; ++iterFeatAlign)
					{
						if (iterFeatAlign->pCDInfo->m_bIsStructDom)
						{
							//dimSDIter[rfidx] = iterFeatAlign;
							//bHasMotifs = true;
							break;
						}
						
						bool bIsSpecific = iterFeatAlign->pAlign->m_bRep && iterFeatAlign->pAlign->m_bSpecQualified;
						const char * pType = ANNOTTYPE_SPECIFIC;
						const list <TOflCDFeat> * pFeats = &(iterFeatAlign->pCDInfo->m_lstSpecFeatures);;
						unsigned int uiSrcPSSMId = iterFeatAlign->pAlign->m_uiPSSMID;
						
						if (!bIsSpecific && NULL != iterFeatAlign->pRootCDInfo)
						{
							pType = ANNOTTYPE_GENERIC;
							pFeats = &(iterFeatAlign->pCDInfo->m_lstGenFeatures);
							uiSrcPSSMId = iterFeatAlign->pCDInfo->m_uiHierarchyRoot;
						}
							
						for (list<TOflCDFeat> :: const_iterator iterFeat = pFeats->begin(); iterFeat != pFeats->end(); ++iterFeat)
						{
							if (TOflCDFeat::eType_StructMotif == iterFeat->m_iType)
							{
								dimMotifTypes[rfidx].push_back(__MotifType(iterFeat, iterFeatAlign, bIsSpecific, uiSrcPSSMId));
								bHasMotifs = true;
								continue;
							}
							
							CSegSet featsegs(*iterFeat);
							iterFeatAlign->pAlign->MapSegSet(featsegs, false);
							
							if (!featsegs.IsEmpty() && ((double)featsegs.GetTotalResidues() / (double)iterFeat->GetCompleteSize() > 0.8))
							{
								vector<CSegSet::TResiduePos> vecMappedPos;
								iterFeatAlign->pAlign->GetTranslatedPosMap(featsegs, vecMappedPos);
      	  	  	
								if (iterFeat->MotifCheck(vecMappedPos, m_dimTranslated[rfidx]) > 0) continue;	//failed motif check
									
								os << idxBlObj << DELIMIT << m_strQueryID << '[' << rf << ']' << DELIMIT << pType << DELIMIT << iterFeat->m_strTitle << DELIMIT;
								
								const TSegs &segs = featsegs.GetSegs();
								
								size_t mapped = 0;
								
								char dimDelimit[2] = {0, 0};
								
								for (TSegs::const_iterator iterSeg = segs.begin(); iterSeg != segs.end(); ++iterSeg)
								{
									int xres = iterFeatAlign->pAlign->Pr2NAPlus(iterSeg->from);
									for (int res = iterSeg->from; res <= iterSeg->to; ++res)
									{
										os << dimDelimit << (char)toupper(m_dimTranslated[rfidx][res]) << (xres + COORDSBASE) << '-';
										xres += RF_SIZE - 1;
										os << (xres + COORDSBASE);
										++xres;
										dimDelimit[0] = COORDELIMIT;
										++mapped;
									}
								}
								os << DELIMIT << iterFeat->GetTotalResidues() << DELIMIT << mapped << DELIMIT << uiSrcPSSMId << endl;
							}
						}
					}
					
					// -- non-motif from sd domains
					for ( ; iterFeatAlignEnd != iterFeatAlign; ++iterFeatAlign)
					{
						//if (iterFeatAlign->pCDInfo->m_bIsStructDom)
						//{
						//	dimSDIter[rfidx] = iterFeatAlign;
						//	bHasMotifs = true;
						//	break;
						//}
						
						//bool bIsSpecific = iterFeatAlign->pAlign->m_bRep && iterFeatAlign->pAlign->m_bSpecQualified;
						//const char * pType = ANNOTTYPE_SPECIFIC;
						//const list <TOflCDFeat> * pFeats = &(iterFeatAlign->pCDInfo->m_lstSpecFeatures);;
						//unsigned int uiSrcPSSMId = iterFeatAlign->pAlign->m_uiPSSMID;
						//
						//if (!bIsSpecific && NULL != iterFeatAlign->pRootCDInfo)
						//{
						//	pType = ANNOTTYPE_GENERIC;
						//	pFeats = &(iterFeatAlign->pCDInfo->m_lstGenFeatures);
						//	uiSrcPSSMId = iterFeatAlign->pCDInfo->m_uiHierarchyRoot;
						//}
							
						for (list<TOflCDFeat> :: const_iterator iterFeat = iterFeatAlign->pCDInfo->m_lstSpecFeatures.begin(); iterFeat != iterFeatAlign->pCDInfo->m_lstSpecFeatures.end(); ++iterFeat)
						{
							if (TOflCDFeat::eType_StructMotif == iterFeat->m_iType)
							{
								dimMotifTypes[rfidx].push_back(__MotifType(iterFeat, iterFeatAlign, true, 0));
								bHasMotifs = true;
								continue;
							}
							
							CSegSet featsegs(*iterFeat);
							iterFeatAlign->pAlign->MapSegSet(featsegs, false);
							
							if (!featsegs.IsEmpty() && ((double)featsegs.GetTotalResidues() / (double)iterFeat->GetCompleteSize() > 0.8))
							{
								vector<CSegSet::TResiduePos> vecMappedPos;
								iterFeatAlign->pAlign->GetTranslatedPosMap(featsegs, vecMappedPos);
      	  	  	
								if (iterFeat->MotifCheck(vecMappedPos, m_dimTranslated[rfidx]) > 0) continue;	//failed motif check
									
								os << idxBlObj << DELIMIT << m_strQueryID << '[' << rf << ']' << DELIMIT << ANNOTTYPE_SPECIFIC << DELIMIT << iterFeat->m_strTitle << DELIMIT;
								
								const TSegs &segs = featsegs.GetSegs();
								
								size_t mapped = 0;
								
								char dimDelimit[2] = {0, 0};
								
								for (TSegs::const_iterator iterSeg = segs.begin(); iterSeg != segs.end(); ++iterSeg)
								{
									int xres = iterFeatAlign->pAlign->Pr2NAPlus(iterSeg->from);
									for (int res = iterSeg->from; res <= iterSeg->to; ++res)
									{
										os << dimDelimit << (char)toupper(m_dimTranslated[rfidx][res]) << (xres + COORDSBASE) << '-';
										xres += RF_SIZE - 1;
										os << (xres + COORDSBASE);
										++xres;
										dimDelimit[0] = COORDELIMIT;
										++mapped;
									}
								}
								os << DELIMIT << iterFeat->GetTotalResidues() << DELIMIT << mapped << DELIMIT << 0/*iterFeatAlign->pAlign->m_uiPSSMID*/ << endl;
							}
						}
					}
				}
				
				for (int rfidx = iNegRF; rfidx < TOTAL_RFS; ++rfidx)	//minus strand reading frames
				{
					int rf = 2 - rfidx;
					vector<TOflCDQuery::__TOflAlignRecord> :: const_iterator iterFeatAlign = dimFeatVecs[rfidx].begin(), iterFeatAlignEnd = dimFeatVecs[rfidx].end();
					for ( ; iterFeatAlignEnd != iterFeatAlign; ++iterFeatAlign)
					{
						if (iterFeatAlign->pCDInfo->m_bIsStructDom)
						{
							//dimSDIter[rfidx] = iterFeatAlign;
							//bHasMotifs = true;
							break;
						}
						
						bool bIsSpecific = iterFeatAlign->pAlign->m_bRep && iterFeatAlign->pAlign->m_bSpecQualified;
						const char * pType = ANNOTTYPE_SPECIFIC;
						const list <TOflCDFeat> * pFeats = &(iterFeatAlign->pCDInfo->m_lstSpecFeatures);;
						unsigned int uiSrcPSSMId = iterFeatAlign->pAlign->m_uiPSSMID;
						
						if (!bIsSpecific && NULL != iterFeatAlign->pRootCDInfo)
						{
							pType = ANNOTTYPE_GENERIC;
							pFeats = &(iterFeatAlign->pCDInfo->m_lstGenFeatures);
							uiSrcPSSMId = iterFeatAlign->pCDInfo->m_uiHierarchyRoot;
						}
							
						for (list<TOflCDFeat> :: const_iterator iterFeat = pFeats->begin(); iterFeat != pFeats->end(); ++iterFeat)
						{
							if (TOflCDFeat::eType_StructMotif == iterFeat->m_iType)
							{
								dimMotifTypes[rfidx].push_back(__MotifType(iterFeat, iterFeatAlign, bIsSpecific, uiSrcPSSMId));
								bHasMotifs = true;
								continue;
							}
							
							CSegSet featsegs(*iterFeat);
							iterFeatAlign->pAlign->MapSegSet(featsegs, false);
							
							if (!featsegs.IsEmpty() && ((double)featsegs.GetTotalResidues() / (double)iterFeat->GetCompleteSize() > 0.8))
							{
								vector<CSegSet::TResiduePos> vecMappedPos;
								iterFeatAlign->pAlign->GetTranslatedPosMap(featsegs, vecMappedPos);
      	  	  	
								if (iterFeat->MotifCheck(vecMappedPos, m_dimTranslated[rfidx]) > 0) continue;	//failed motif check
									
								os << idxBlObj << DELIMIT << m_strQueryID << '[' << rf << ']' << DELIMIT << pType << DELIMIT << iterFeat->m_strTitle << DELIMIT;
								
								const TSegs &segs = featsegs.GetSegs();
								
								size_t mapped = 0;
								
								char dimDelimit[2] = {0, 0};
								
								for (TSegs::const_iterator iterSeg = segs.begin(); iterSeg != segs.end(); ++iterSeg)
								{
									int xres = iterFeatAlign->pAlign->Pr2NAMinus(iterSeg->from);
									for (int res = iterSeg->from; res <= iterSeg->to; ++res)
									{
										os << dimDelimit << (char)toupper(m_dimTranslated[rfidx][res]) << (xres + COORDSBASE) << '-';
										xres -= RF_SIZE - 1;
										os << (xres + COORDSBASE);
										--xres;
										dimDelimit[0] = COORDELIMIT;
										++mapped;
									}
								}
								os << DELIMIT << iterFeat->GetTotalResidues() << DELIMIT << mapped << DELIMIT << uiSrcPSSMId << endl;
							}
						}
					}
					
					// -- non-motif features from SD domain
					for ( ; iterFeatAlignEnd != iterFeatAlign; ++iterFeatAlign)
					{
						//if (iterFeatAlign->pCDInfo->m_bIsStructDom)
						//{
						//	//dimSDIter[rfidx] = iterFeatAlign;
						//	//bHasMotifs = true;
						//	break;
						//}
						
						//bool bIsSpecific = iterFeatAlign->pAlign->m_bRep && iterFeatAlign->pAlign->m_bSpecQualified;
						//const char * pType = ANNOTTYPE_SPECIFIC;
						//const list <TOflCDFeat> * pFeats = &(iterFeatAlign->pCDInfo->m_lstSpecFeatures);;
						//unsigned int uiSrcPSSMId = iterFeatAlign->pAlign->m_uiPSSMID;
						//
						//if (!bIsSpecific && NULL != iterFeatAlign->pRootCDInfo)
						//{
						//	pType = ANNOTTYPE_GENERIC;
						//	pFeats = &(iterFeatAlign->pCDInfo->m_lstGenFeatures);
						//	uiSrcPSSMId = iterFeatAlign->pCDInfo->m_uiHierarchyRoot;
						//}
							
						for (list<TOflCDFeat> :: const_iterator iterFeat = iterFeatAlign->pCDInfo->m_lstSpecFeatures.begin(); iterFeat != iterFeatAlign->pCDInfo->m_lstSpecFeatures.end(); ++iterFeat)
						{
							if (TOflCDFeat::eType_StructMotif == iterFeat->m_iType)
							{
								dimMotifTypes[rfidx].push_back(__MotifType(iterFeat, iterFeatAlign, true, 0));
								bHasMotifs = true;
								continue;
							}
							
							CSegSet featsegs(*iterFeat);
							iterFeatAlign->pAlign->MapSegSet(featsegs, false);
							
							if (!featsegs.IsEmpty() && ((double)featsegs.GetTotalResidues() / (double)iterFeat->GetCompleteSize() > 0.8))
							{
								vector<CSegSet::TResiduePos> vecMappedPos;
								iterFeatAlign->pAlign->GetTranslatedPosMap(featsegs, vecMappedPos);
      	  	  	
								if (iterFeat->MotifCheck(vecMappedPos, m_dimTranslated[rfidx]) > 0) continue;	//failed motif check
									
								os << idxBlObj << DELIMIT << m_strQueryID << '[' << rf << ']' << DELIMIT << ANNOTTYPE_SPECIFIC << DELIMIT << iterFeat->m_strTitle << DELIMIT;
								
								const TSegs &segs = featsegs.GetSegs();
								
								size_t mapped = 0;
								
								char dimDelimit[2] = {0, 0};
								
								for (TSegs::const_iterator iterSeg = segs.begin(); iterSeg != segs.end(); ++iterSeg)
								{
									int xres = iterFeatAlign->pAlign->Pr2NAMinus(iterSeg->from);
									for (int res = iterSeg->from; res <= iterSeg->to; ++res)
									{
										os << dimDelimit << (char)toupper(m_dimTranslated[rfidx][res]) << (xres + COORDSBASE) << '-';
										xres -= RF_SIZE - 1;
										os << (xres + COORDSBASE);
										--xres;
										dimDelimit[0] = COORDELIMIT;
										++mapped;
									}
								}
								os << DELIMIT << iterFeat->GetTotalResidues() << DELIMIT << mapped << DELIMIT << 0/*uiSrcPSSMId*/ << endl;
							}
						}
					}
				}	//negative strand reading frames
				
				if (bHasRegFeats)
					os << FEATEND << endl;
					
				if (bHasMotifs)
				{
					os << MOTIFSTART << endl;
					
					for (int rfidx = 0; rfidx < iNegRF; ++rfidx)
					{
						int rf = rfidx + 1;
						for (vector<__MotifType> :: const_iterator iterM = dimMotifTypes[rfidx].begin(), iterMEnd = dimMotifTypes[rfidx].end(); iterMEnd != iterM; ++iterM)
						{
							CSegSet featsegs(*iterM->iterMotifFeat);
							iterM->iterAlignRec->pAlign->MapSegSet(featsegs);
							
							if (!featsegs.IsEmpty() && ((double)featsegs.GetTotalResidues() / (double)iterM->iterMotifFeat->GetCompleteSize() > 0.8))
							{
								os << idxBlObj << DELIMIT << m_strQueryID << '[' << rf << ']' << DELIMIT << iterM->iterMotifFeat->m_strTitle << DELIMIT;
								int res0 = featsegs.GetLeft(), res1 = featsegs.GetRight();
								os << res0 + COORDSBASE << DELIMIT << res1 + RF_SIZE - 1 + COORDSBASE << DELIMIT << iterM->uiSrcPSSMId << endl;
							}
						}
						
						// -- now start SD
						//vector<TOflCDQuery::__TOflAlignRecord> :: const_iterator iterFeatAlign = dimSDIter[rfidx], iterFeatAlignEnd = dimFeatVecs[rfidx].end();
						//for ( ; iterFeatAlignEnd != iterFeatAlign; ++iterFeatAlign)
						//{
						//	//int rfidx = iterFeatAlign->pAlign->GetRFIdx();
						//	for (list<TOflCDFeat> :: const_iterator iterFeat = iterFeatAlign->pCDInfo->m_lstSpecFeatures.begin(); iterFeat != iterFeatAlign->pCDInfo->m_lstSpecFeatures.end(); ++iterFeat)
						//	{
						//		CSegSet featsegs(*iterFeat);
						//		iterFeatAlign->pAlign->MapSegSet(featsegs);
						//		if (!featsegs.IsEmpty())
						//		{
						//			os << idxBlObj << DELIMIT << m_strQueryID << '[' << rf << ']' << DELIMIT << iterFeat->m_strTitle << DELIMIT;
						//			//int rfidx = iterFeatAlign->pAlign->GetRFIdx();
						//			int res0 = featsegs.GetLeft(), res1 = featsegs.GetRight();
						//			os << res0 + COORDSBASE << DELIMIT << res1 + RF_SIZE - 1 + COORDSBASE << DELIMIT << '-' << endl;
						//		}
						//	}
						//}
					}
					
					// -- minus strand
					for (int rfidx = iNegRF; rfidx < TOTAL_RFS; ++rfidx)
					{
						int rf = 2 - rfidx;
						for (vector<__MotifType> :: const_iterator iterM = dimMotifTypes[rfidx].begin(), iterMEnd = dimMotifTypes[rfidx].end(); iterMEnd != iterM; ++iterM)
						{
							CSegSet featsegs(*iterM->iterMotifFeat);
							iterM->iterAlignRec->pAlign->MapSegSet(featsegs);
							
							if (!featsegs.IsEmpty() && ((double)featsegs.GetTotalResidues() / (double)iterM->iterMotifFeat->GetCompleteSize() > 0.8))
							{
								os << idxBlObj << DELIMIT << m_strQueryID << '[' << rf << ']' << DELIMIT << iterM->iterMotifFeat->m_strTitle << DELIMIT;
								int res0 = featsegs.GetLeft(), res1 = featsegs.GetRight();
								os << res1 + COORDSBASE << DELIMIT << res0 - RF_SIZE + 1 + COORDSBASE << DELIMIT << iterM->uiSrcPSSMId << endl;
							}
						}
						
						// -- now start SD
						//vector<TOflCDQuery::__TOflAlignRecord> :: const_iterator iterFeatAlign = dimSDIter[rfidx], iterFeatAlignEnd = dimFeatVecs[rfidx].end();
						//for ( ; iterFeatAlignEnd != iterFeatAlign; ++iterFeatAlign)
						//{
						//	//int rfidx = iterFeatAlign->pAlign->GetRFIdx();
						//	for (list<TOflCDFeat> :: const_iterator iterFeat = iterFeatAlign->pCDInfo->m_lstSpecFeatures.begin(); iterFeat != iterFeatAlign->pCDInfo->m_lstSpecFeatures.end(); ++iterFeat)
						//	{
						//		CSegSet featsegs(*iterFeat);
						//		iterFeatAlign->pAlign->MapSegSet(featsegs);
						//		if (!featsegs.IsEmpty())
						//		{
						//			os << idxBlObj << DELIMIT << m_strQueryID << '[' << rf << ']' << DELIMIT << iterFeat->m_strTitle << DELIMIT;
						//			//int rfidx = iterFeatAlign->pAlign->GetRFIdx();
						//			int res0 = featsegs.GetLeft(), res1 = featsegs.GetRight();
						//			os << os << res1 + COORDSBASE << DELIMIT << res0 - RF_SIZE + 1 + COORDSBASE << DELIMIT << '-' << endl;
						//		}
						//	}
						//}
					}
					os << MOTIFEND << endl;
				}
			}
		}
	}
	
	os << QUERYEND << endl;
	if (!g_bSilent) cerr << "End processing Query " << m_strDefline << endl;
}


/**********************************************************************
*	Biodata structure -- The main processor
***********************************************************************/
class COflRpsbPostProcessor: public COflDomClstInfo
{
public:
	COflRpsbPostProcessor(void): COflDomClstInfo() {};
	
	// -- return count of objects processed
	int ProcessRpsbDataStream(istream &is, ostream &os, double evcut, int mode) const;
	
private:
	
	struct TCDSortFacility
	{
		TOflCDAlignInfo *pAlign;
		size_t ulIdx;
		
		// -- added 12/03/2013 for sorting modification
		const TOflCDInfo *pCDInfo;
		size_t ulGapsIdx;
		
		
		TCDSortFacility(void): pAlign(nullptr), ulIdx(0), pCDInfo(nullptr), ulGapsIdx(string::npos) {};
	};
	
	struct TSortCDByEValue
	{
		bool operator () (const TCDSortFacility &p1, const TCDSortFacility &p2)
		{
			if (p1.pAlign->m_dEValue < p2.pAlign->m_dEValue) return true;
			else if (p1.pAlign->m_dEValue > p2.pAlign->m_dEValue) return false;
			else return (p1.pAlign->m_dBitScore > p2.pAlign->m_dBitScore);
		}
	};
	
	void x_Calculate(vector<TOflCDAlignInfo> &aligns, const vector<size_t> &selIdx, TOflAlignIndice &dst) const;
	static void x_SortReadingFrames(vector<size_t> rfIndice[TOTAL_RFS], TOflCdQueryEx &rTarget);
	static const char * const DOMSRCSIGS[];
	static const int DOMSRCCNTS[];
	static int DomAccType(const string &acc);
};

const char * const COflRpsbPostProcessor::DOMSRCSIGS[] = {"CD", "CHL", "COG", "MTH", "PFAM", "PHA", "PLN", "PRK", "PTZ", "SMART", "TIGR", NULL};
const int COflRpsbPostProcessor::DOMSRCCNTS[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0};

int COflRpsbPostProcessor::DomAccType(const string &acc)
{
	int iSig = 0, iChar = 0, iAccChar = 0;
	
	int acclen = (int)acc.size();
	while (nullptr != DOMSRCSIGS[iSig] && 0 != DOMSRCSIGS[iSig][iChar] && iAccChar < acclen)
	{
		char accChar = acc[iAccChar];
		if ('a' <= accChar && 'z' >= accChar) accChar -= 0x20;	//to uppercase
		if (accChar == DOMSRCSIGS[iSig][iChar])
		{
			++iAccChar;
			++iChar;
		}
		else
		{

			++iSig;
			iChar = 0;
			iAccChar = 0;
		}
	}
	return iSig;
}


int COflRpsbPostProcessor::ProcessRpsbDataStream(istream &is, ostream &os, double evcut, int mode) const
{
	int objCounter = 0;
	while (is.good())
	{
		CBlastOutput objBlastOutput;
		try
		{
			ObjStreamIn<CBlastOutput> (is, objBlastOutput, eSerial_Xml);
		}
		catch (...)
		{
			return objCounter;
		}
		// -- start to output
		const CBlastOutput::TParam &params = objBlastOutput.GetParam();
		os << SESSIONSTART << DELIMIT << ++objCounter << DELIMIT << objBlastOutput.GetVersion() << DELIMIT << objBlastOutput.GetDb() << DELIMIT << params.GetMatrix() << DELIMIT << params.GetExpect() << endl;
		if (!g_bSilent) cerr << "RPSBlast session " << objCounter << " start..." << endl;
		const CBlastOutput::TIterations &rIters = objBlastOutput.GetIterations();

		for (CBlastOutput::TIterations::const_iterator iter = rIters.begin(), iterEnd = rIters.end(); iterEnd != iter; ++iter)
		{
			TOflCdQueryEx cdqVal;
			cdqVal.ProcessQueryResult(*iter, evcut);
			if (!g_bSilent) cerr << "Calculating...." << endl;
			if (cdqVal.m_bIsProtein)
			{
				size_t ttlAligns = cdqVal.m_vecAlignments.size();


				vector<size_t> ttlIndice;
				ttlIndice.reserve(ttlAligns);
				for (size_t i = 0; i < ttlAligns; ++i)
					ttlIndice.push_back(i);
				x_Calculate(cdqVal.m_vecAlignments, ttlIndice, cdqVal);
			}
			else
			{
				size_t ttlAligns = cdqVal.m_vecAlignments.size();
				vector<size_t> rfIndice[TOTAL_RFS];
				for (size_t i = 0; i < TOTAL_RFS; ++i)
					rfIndice[i].reserve(ttlAligns);
					
				x_SortReadingFrames(rfIndice, cdqVal);
				
				for (size_t i = 0; i < TOTAL_RFS; ++i)
					x_Calculate(cdqVal.m_vecAlignments, rfIndice[i], cdqVal.m_dimSplitAligns[i]);
			}
			
			cdqVal.Print(objCounter, *this, os, mode);
		}
		os << SESSIONEND << DELIMIT << objCounter << endl;
		if (!g_bSilent) cerr << "RPSBlast session " << objCounter << " done..." << endl;
	}
	return objCounter;
}

void COflRpsbPostProcessor::x_SortReadingFrames(vector<size_t> rfIndice[TOTAL_RFS], TOflCdQueryEx &rTarget)
{
	for (size_t idx = 0, len = rTarget.m_vecAlignments.size(); idx < len; ++idx)
		rfIndice[rTarget.m_vecAlignments[idx].GetRFIdx()].push_back(idx);
}

void COflRpsbPostProcessor::x_Calculate(vector<TOflCDAlignInfo> &aligns, const vector<size_t> &selIdx, TOflAlignIndice &dst) const
//This implementation does not promote NCBI-curated models
{
	vector <CSegSet> vecGaps;
	vecGaps.push_back(CSegSet());	//sentinel
	size_t gapIdx = 0;
	vector <TCDSortFacility> vecPrivileged, vecNonMulti, vecLongMultiDom, vecMultiDom, vecConcise, vecSDs;

	if (!selIdx.empty())
	{
		TCDSortFacility dummy;
		TSortCDByEValue sortingEV;
		
		for (vector<size_t> :: const_iterator iterIdx = selIdx.begin(), iterIdxEnd = selIdx.end(); iterIdx != iterIdxEnd; ++iterIdx)
		{
			TOflCDAlignInfo &rAlign = aligns[*iterIdx];
			const TOflCDInfo *pCDInfo = FindCDInfo(rAlign.m_uiPSSMID);
			
			if (NULL == pCDInfo)	//unrecognized CD
			{
				cerr << "PSSMId " << rAlign.m_uiPSSMID << " not found. Ignored" << endl;
				continue;
			}
			rAlign.m_dSeqIdentity = (double)(rAlign.m_iNumIdent) / (double)(pCDInfo->m_uiLength) * 100.0;
			
			int iNMissing = *(rAlign.m_vecSStarts.begin());
			int iCMissing = pCDInfo->m_uiLength - (*(rAlign.m_vecSStarts.rbegin()) + *(rAlign.m_vecLens.rbegin()) - 1);
			
			rAlign.m_dNMissing = (double)(iNMissing) / (double)(pCDInfo->m_uiLength);
			rAlign.m_dCMissing = (double)(iCMissing) / (double)(pCDInfo->m_uiLength);
			rAlign.m_uiAlignedLen = pCDInfo->m_uiLength - iNMissing - iCMissing;
			rAlign.m_dAlignedPct = (double)(rAlign.m_uiAlignedLen) / (double)(pCDInfo->m_uiLength) * 100.0;
			rAlign.m_bSpecQualified = (pCDInfo->m_dMinBitScore > 0.0 ? rAlign.m_dBitScore >= pCDInfo->m_dMinBitScore : !pCDInfo->m_bCurated);
			// -- sorting
			
			dummy.pAlign = &(rAlign);
			dummy.ulIdx = *iterIdx;
			dummy.pCDInfo = pCDInfo;
			dummy.ulGapsIdx = string::npos;
			
			rAlign.CalcMasterGaps(TOflCDAlignInfo::GAP_THRESHOLD, vecGaps[gapIdx]);
			if (!vecGaps[gapIdx].IsEmpty())
			{
				dummy.ulGapsIdx = gapIdx;
				++gapIdx;
				vecGaps.push_back(CSegSet());
			}
		
			if (pCDInfo->m_bIsStructDom)
			{
				rAlign.m_iRepClass = 0x2;
				vecSDs.push_back(dummy);
			}
			else if (pCDInfo->m_bCurated || !pCDInfo->m_bMultiDom)	//non-curated single dom
			{
				rAlign.m_iRepClass = 0;
				vecNonMulti.push_back(dummy);
			}
			//else if (rAlign.m_dAlignedPct > 90.0)	//long
			//{
			//	rAlign.m_iRepClass = 1;
			//	vecLongMultiDom.push_back(dummy);
			//}
			else
			{
				rAlign.m_iRepClass = 1;
				vecMultiDom.push_back(dummy);
			}
		}
	  map<int, TDomSrcCount> dimAccTypeCount, dimAccTypeCount2;
		
		// -- sort
		sort(vecNonMulti.begin(), vecNonMulti.end(), sortingEV);
		for (vector <TCDSortFacility>::const_iterator iter = vecNonMulti.begin(), iterEnd = vecNonMulti.end(); iterEnd != iter; ++iter)
		{

			CSegSet s_overlaps;	// to calculate combined overlap region
			
			CSegSet thisSegs;
			thisSegs.AddSeg(iter->pAlign->m_iFrom, iter->pAlign->m_iTo);
			if (string::npos != iter->ulGapsIdx)
				thisSegs.Clip(vecGaps[iter->ulGapsIdx]);
			
			int iThisLength = (int)thisSegs.GetTotalResidues();
			double dThisLength = (double)iThisLength;
			
			for (vector <TCDSortFacility>::const_iterator iterRep = vecConcise.begin(), iterRepEnd = vecConcise.end(); iterRepEnd != iterRep; ++iterRep)
			{
				CSegSet repSegs;
				repSegs.AddSeg(iterRep->pAlign->m_iFrom, iterRep->pAlign->m_iTo);
				
				if (string::npos != iterRep->ulGapsIdx)
					repSegs.Clip(vecGaps[iterRep->ulGapsIdx]);
				
				int iRepLength = (int)repSegs.GetTotalResidues();
				double dRepLength = (double)iRepLength;
				
				// -- find gaps. any gap >= half of the domain model length are considered a gap and excluded from overlapping
				
				CSegSet olSegs(thisSegs);
				olSegs.Cross(repSegs);
				
				int iOverlapLen = (int)olSegs.GetTotalResidues();
				double dOverlapLen = (double)iOverlapLen;

				if (dOverlapLen > 0)
				{
				
					if (dOverlapLen / dThisLength > 0.5 || (dOverlapLen / dRepLength > 0.5 && iter->pCDInfo->m_uiClusterPSSMID == iterRep->pCDInfo->m_uiClusterPSSMID))	//mutually overlap > 0.5
					//if (dOverlapLen / dThisLength + dOverlapLen / dRepLength > 1.0)	//mutually overlap > 0.5
					{
						iter->pAlign->m_iRegionIdx = iterRep->pAlign->m_iRegionIdx;
						iter->pAlign->m_bRep = false;
						
						goto labelNextNonMulti;
					}
				}
				s_overlaps.Merge(olSegs);

			}
			
			if ((iter->pAlign->m_bRep = (((double)(s_overlaps.GetTotalResidues()) / dThisLength) < 0.5)))	//new region
			{
				iter->pAlign->m_iRegionIdx = vecConcise.size();
				vecConcise.push_back(*iter);

			}
		labelNextNonMulti:
			
			map<int, TDomSrcCount> :: iterator iterSrcCounter = dimAccTypeCount.emplace(iter->pAlign->m_iRegionIdx, TDomSrcCount()).first;
			if (iterSrcCounter->second.CountSrc(iter->pCDInfo->m_strAccession))
			{
				if (iter->pCDInfo->m_bCurated)
				{
					dst.m_vecQualifiedFeatIndice.push_back(iter->ulIdx);
				}
				dst.m_vecStdIndice.push_back(iter->ulIdx);
				iter->pAlign->m_bRep = true;
			}
			dst.m_vecSortedIndice.push_back(iter->ulIdx);
		}
		
		// -- When done, convert setup dst.m_vecConciseIndice
		size_t ttlConcise = vecConcise.size();
		if (ttlConcise > 0)
		{
			dst.m_vecConciseIndice.reserve(ttlConcise);
			for (size_t i = 0; i < ttlConcise; ++i)
				dst.m_vecConciseIndice.push_back(vecConcise[i].ulIdx);
		}
		// -- now deal with multi. first process long (> 90% multi)
		//sort(vecLongMultiDom.begin(), vecLongMultiDom.end(), sortingEV);
		//vecPrivileged.clear();
		//for (vector <TCDSortFacility>::const_iterator iter = vecLongMultiDom.begin(), iterEnd = vecLongMultiDom.end(); iterEnd != iter; ++iter)
		//{
		//	CSegSet thisSegs;
		//	thisSegs.AddSeg(iter->pAlign->m_iFrom, iter->pAlign->m_iTo);
		//	if (string::npos != iter->ulGapsIdx)
		//		thisSegs.Clip(vecGaps[iter->ulGapsIdx]);
		//	
		//	int iThisLength = (int)thisSegs.GetTotalResidues();
		//	double dThisLength = (double)iThisLength;
		//	
		//	
		//	CSegSet m_overlaps;	//only consider multi overlaping
		//	double dMOverlap = 0.0;
		//	
		//	for (vector <TCDSortFacility>::const_iterator iterRep = vecPrivileged.begin(), iterRepEnd = vecPrivileged.end(); iterRepEnd != iterRep; ++iterRep)
		//	{
		//		if (!(iterRep->pAlign->m_iFrom > iter->pAlign->m_iTo || iterRep->pAlign->m_iTo < iter->pAlign->m_iFrom))	//overlap
		//		{
		//			CSegSet repSegs;
		//			repSegs.AddSeg(iterRep->pAlign->m_iFrom, iterRep->pAlign->m_iTo);
		//		
		//			if (string::npos != iterRep->ulGapsIdx)
		//				repSegs.Clip(vecGaps[iterRep->ulGapsIdx]);
		//		
		//			int iRepLength = (int)repSegs.GetTotalResidues();
		//			double dRepLength = (double)iRepLength;
		//			
		//			
		//			CSegSet olSegs(thisSegs);
		//			olSegs.Cross(repSegs);
		//			
		//			int iOverlapLen = (int)olSegs.GetTotalResidues();
		//			
		//			double dOverlapLen = (double)(iOverlapLen);
		//			
		//			if (dOverlapLen / dThisLength + dOverlapLen / dRepLength > 1.0)	//mutually overlap > 0.5
		//			{
		//				iter->pAlign->m_iRegionIdx = iterRep->pAlign->m_iRegionIdx;
		//				iter->pAlign->m_bRep = false;
		//				vecMultiDom.push_back(*iter);
		//				goto labelNextLongMulti;
		//			}
		//			m_overlaps.Merge(olSegs);
		//		}
		//	}
		//	
		//	dMOverlap = (double)(m_overlaps.GetTotalResidues()) / dThisLength;
		//	
		//	if (dMOverlap <= 0.5)	//new multi-dom region -- rep!
		//	{
		//		iter->pAlign->m_iRegionIdx = dst.m_vecConciseIndice.size();
		//		iter->pAlign->m_bRep = true;
		//		
		//		map<int, TDomSrcCount> :: iterator iterSrcCounter = dimAccTypeCount2.emplace(iter->pAlign->m_iRegionIdx, TDomSrcCount()).first;
		//		iterSrcCounter->second.CountSrc(iter->pCDInfo->m_strAccession);
		//		dst.m_vecConciseIndice.push_back(iter->ulIdx);
		//		vecConcise.push_back(*iter);
		//		dst.m_vecStdIndice.push_back(iter->ulIdx);
		//		dst.m_vecSortedIndice.push_back(iter->ulIdx);
		//		vecPrivileged.push_back(*iter);
		//	}
		//	else vecMultiDom.push_back(*iter);
		//		
		//labelNextLongMulti:;
		//}
		// -- now deal with short multi
		// -- if no concise at all, then vecPrivileged should be already empty. so it will not change the status that vecPrivileged always contains the "concise" of multidoms to display if no concise found.
		
		vecPrivileged.clear();
		sort(vecMultiDom.begin(), vecMultiDom.end(), sortingEV);
		for (vector <TCDSortFacility>::const_iterator iter = vecMultiDom.begin(), iterEnd = vecMultiDom.end(); iterEnd != iter; ++iter)
		{
			CSegSet s_overlaps;	// to calculate combined overlap region
			CSegSet m_overlaps;
			
			double dSOverlap = 0.0;	//put here for scope reason
			double dMOverlap = 0.0;
			
			CSegSet thisSegs;
			thisSegs.AddSeg(iter->pAlign->m_iFrom, iter->pAlign->m_iTo);
			if (string::npos != iter->ulGapsIdx)
				thisSegs.Clip(vecGaps[iter->ulGapsIdx]);
			
			int iThisLength = (int)thisSegs.GetTotalResidues();
			double dThisLength = (double)iThisLength;
			
			//double dThisLength = (double)(iter->pAlign->m_iTo - iter->pAlign->m_iFrom + 1);
			for (vector <TCDSortFacility>::const_iterator iterRep = vecConcise.begin(), iterRepEnd = vecConcise.end(); iterRepEnd != iterRep; ++iterRep)
			{
				if (!(iterRep->pAlign->m_iFrom > iter->pAlign->m_iTo || iterRep->pAlign->m_iTo < iter->pAlign->m_iFrom))	//overlap
				{
					CSegSet repSegs;
					repSegs.AddSeg(iterRep->pAlign->m_iFrom, iterRep->pAlign->m_iTo);
					
					if (string::npos != iterRep->ulGapsIdx)
						repSegs.Clip(vecGaps[iterRep->ulGapsIdx]);
					
					int iRepLength = (int)repSegs.GetTotalResidues();
					double dRepLength = (double)iRepLength;
					
					// -- find gaps. any gap >= half of the domain model length are considered a gap and excluded from overlapping
					
					CSegSet olSegs(thisSegs);
					olSegs.Cross(repSegs);
					
					int iOverlapLen = (int)olSegs.GetTotalResidues();
					
					double dOverlapLen = (double)iOverlapLen;
					
					if (0 == iterRep->pAlign->m_iRepClass)	//single
					{
						s_overlaps.Merge(olSegs);
					}
					else	//multi
					{
						if (dOverlapLen / dThisLength + dOverlapLen / dRepLength > 1.0)	//mutually overlap > 0.5
						{
							iter->pAlign->m_iRegionIdx = iterRep->pAlign->m_iRegionIdx;
							iter->pAlign->m_bRep = false;
							goto labelNextMulti;
						}
						m_overlaps.Merge(olSegs);
					}
				}
			}
			
			
			for (vector <TCDSortFacility>::const_iterator iterRep = vecPrivileged.begin(), iterRepEnd = vecPrivileged.end(); iterRepEnd != iterRep; ++iterRep)
			{
				if (!(iterRep->pAlign->m_iFrom > iter->pAlign->m_iTo || iterRep->pAlign->m_iTo < iter->pAlign->m_iFrom))	//overlap
				{
					
					CSegSet repSegs;
					repSegs.AddSeg(iterRep->pAlign->m_iFrom, iterRep->pAlign->m_iTo);
					
					if (string::npos != iterRep->ulGapsIdx)
						repSegs.Clip(vecGaps[iterRep->ulGapsIdx]);
					
					int iRepLength = (int)repSegs.GetTotalResidues();
					double dRepLength = (double)iRepLength;
					
					// -- find gaps. any gap >= half of the domain model length are considered a gap and excluded from overlapping
					
					CSegSet olSegs(thisSegs);
					olSegs.Cross(repSegs);
					
					int iOverlapLen = (int)olSegs.GetTotalResidues();
					
					double dOverlapLen = (double)iOverlapLen;
					
					if (dOverlapLen / dThisLength + dOverlapLen / dRepLength > 1.0)	//mutually overlap > 0.5
					{
						iter->pAlign->m_iRegionIdx = iterRep->pAlign->m_iRegionIdx;
						iter->pAlign->m_bRep = false;
						goto labelNextMulti;
					}
					m_overlaps.Merge(olSegs);
				}
			}
			
			dSOverlap = (double)(s_overlaps.GetTotalResidues()) / dThisLength;
			dMOverlap = (double)(m_overlaps.GetTotalResidues()) / dThisLength;

			if (dMOverlap <= 0.5)	//new multi-dom region
			{
				if ((dSOverlap <= 0.5 && iter->pAlign->m_dAlignedPct > 50.0) || iter->pAlign->m_bSpecQualified)	//rep!
				{
					iter->pAlign->m_iRegionIdx = dst.m_vecConciseIndice.size() + vecPrivileged.size();
					iter->pAlign->m_bRep = true;
					
					
					dst.m_vecConciseIndice.push_back(iter->ulIdx);
					vecConcise.push_back(*iter);
				}
				vecPrivileged.push_back(*iter);
			}
		labelNextMulti:
			
			map<int, TDomSrcCount> :: iterator iterSrcCounter = dimAccTypeCount2.emplace(iter->pAlign->m_iRegionIdx, TDomSrcCount()).first;
			if (iterSrcCounter->second.CountSrc(iter->pCDInfo->m_strAccession))
				dst.m_vecStdIndice.push_back(iter->ulIdx);
			
			dst.m_vecSortedIndice.push_back(iter->ulIdx);
		
		}
		// -- to avoid empty concise
		
		if (!dst.m_vecSortedIndice.empty())
		{
			if (dst.m_vecConciseIndice.empty())
			{
				size_t iEnd = vecPrivileged.size();
				dst.m_vecConciseIndice.reserve(iEnd);
				for (size_t i = 0; i < iEnd; ++i)
				{
					dst.m_vecConciseIndice.push_back(vecPrivileged[i].ulIdx);
					vecPrivileged[i].pAlign->m_bRep = true;
				}
			}
			
			if (dst.m_vecStdIndice.empty())
			{
				dst.m_vecStdIndice = dst.m_vecSortedIndice;
			}
		}



		// -- finally calculate SDs
		// -- SDs: models are sorted according to evalues, each model reserve their regions on the query sequence. that would need to change the from/to and the aligned segments. Features are 
		// -- then mapped to the query from these regions -- to guarantee non-redundency. Proposed by Aron.
		// -- implementation: calculate a restriction segset for SD models, which will be used later to trim the sites first
		if (!vecSDs.empty())
		{
			sort(vecSDs.begin(), vecSDs.end(), sortingEV);
			//vecConcise.clear();
			dst.m_vecSDIndice.clear();
			dst.m_vecSDIndice.reserve(vecSDs.size());
			CSegSet covered_regions;
			
			for (vector <TCDSortFacility>::const_iterator iter = vecSDs.begin(), iterEnd = vecSDs.end(); iterEnd != iter; ++iter)
			{

				CSegSet sdSegSet;
				sdSegSet.AddSeg(0, iter->pCDInfo->m_uiLength - 1);


				iter->pAlign->MapSegSet(sdSegSet);
				
				sdSegSet.Clip(covered_regions);
				
				if (!sdSegSet.IsEmpty())
				{
					iter->pAlign->m_ClipSet.Clear();
					const TSegs& sdSegs = sdSegSet.GetSegs();

					
					for (TSegs::const_iterator iterSeg = sdSegs.begin(), iterSegEnd = sdSegs.end(); iterSegEnd != iterSeg; ++iterSeg)
					{
						iter->pAlign->m_ClipSet.AddSeg(iterSeg->ori_from, sdSegSet.GetOriTo(iterSeg));

					}
					
					covered_regions.Merge(sdSegSet);
					dst.m_vecSDIndice.push_back(iter->ulIdx);
				}
			}
		}
	}
}

/**********************************************************************
*	Main functions
***********************************************************************/

void PrintUsage(const char * pname, ostream &os)
{
	os << PROGRAM_TITLE << endl;
	os << "This utility processes domain hit data produced by local RPS-BLAST and" << endl;
	os << "and generates domain family and/or superfamily annotations on the query sequences." << endl;
	os << "Some data files (all downloadable from NCBI ftp site) are required for this program" << endl;
	os << "to function correctly. " << endl << endl;
	os << "Usage:" << endl;
	os << pname << " [-c<config-file>] [-i<in-filename>] [-o<out-filename>] [-e<evalue-cutoff>] [-m<data-mode>] [-d<target-data>] [-h]" << endl << endl;
	os << "\t-c<config-file>" << endl;
	os << "\t\tBy default, the program looks for a <progname>.ini file in current directory" << endl;
	os << "\t\tand loads configuration data from it. This -c switch can be used to override the default and load" << endl;
	os << "\t\tconfiguration data from a different file" << endl << endl;
	os << "\t-i<in-filename>" << endl;
	os << "\t\tThe file contains data generated by local RPS-BLAST program. If omitted, default to stdin" << endl << endl;
	os << "\t-o<out-filename>" << endl;
	os << "\t\tFile to which the processed results should be written. If omitted, default to stdout" << endl << endl;
	os << "\t-e<evalue-cutoff>" << endl;
	os << "\t\tOnly processes hits with evalues better (smaller in value) than the designated." << endl;
	os << "\t\tIf omitted, default to 0.01" << endl << endl;
	os << "\t-m<data-mode>" << endl;
	os << "\t\tSelect redundancy level of domain hit data. Valid options are \"rep\" (concise), \"std\"(standard)" << endl;
	os << "\t\tand \"full\" (all hits)" << endl << endl;
	os << "\t-t<target-data>" << endl;
	os << "\t\tSelect desired (target) data. Valid options are \"doms\", \"feats\" or \"both\". . If omitted, default to \"both\"" << endl << endl;
	os << "\t-d<datapath>" << endl;
	os << "\t\tLocation of data files. Usually data files are set in the [datapath] secion of the configure file (" << pname << ".ini)" << endl;
	os << "\t\tand have the freedom to use any file names. This switch overrides thse settings and looks for the exact file set under the specified path:" << endl;
	os << endl;
	os << "\t\t\t" << CDIDFILE << endl;
	os << "\t\t\t" << CDTRACKFILE << endl;
	os << "\t\t\t" << CLSTLINKFILE << endl;
	os << "\t\t\t" << SPFEATFILE << endl;
	os << "\t\t\t" << GENFEATFILE << endl;
	os << "\t\t\t" << MINBSCOREFILE << endl;
	os << endl;
	os << "\t\tSo be sure to check the data file names after downloading from our ftp site and rename them if this switch is to be used." << endl << endl;
	os << "\t-f" << endl;
	os << "\t\tShow corresponding superfamily information for domain hits. A " << FAMSTART << '/' << FAMEND << " section will appear after" << endl;
	os << "\t\tthe " << DOMSTART << '/' << DOMEND << " section to list superfamily hits corresponding to the representative domain hits in each region." << endl << endl;
	os << "\t-q" << endl;
	os << "\t\tQuiet mode -- do not display program information and progress on stderr" << endl << endl;
	os << "\t-h" << endl;
	os << "\t\tDisplay this help screen, ignore other switches" << endl << endl;
	
	os << "Configuration file:" << endl;
	os << "The [datapath] section contains paths to required data files. All data files must be downloaded from NCBI FTP site:" << endl << endl;
	os << "\tftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/" << endl << endl;
	os << "to the local file system and decompressed. Edit the configuration file so that the paths point to the" << endl;
	os << "correct location of these files" << endl << endl;
	os << "Thanks for using our services and programs." << endl;
}

void ProcessCmdline(int sw, const char *arg)
{

	switch (sw)
	{
	case TRpsbProcCmds::eCfgFile:
		if (NULL == arg)
			cerr << "missing parameter for switch -" << TRpsbProcCmds::m_dimRpsbProcSwitches[TRpsbProcCmds::eCfgFile].m_cToken << ", ignored.." << endl;
		else
			g_strCfgFile = arg;
		break;
	case TRpsbProcCmds::eInFile:
		if (NULL == arg)
			cerr << "missing parameter for switch -" << TRpsbProcCmds::m_dimRpsbProcSwitches[TRpsbProcCmds::eInFile].m_cToken << ", ignored.." << endl;
		else
			g_strSrcFile = arg;
		break;
	case TRpsbProcCmds::eOutFile:
		if (NULL == arg)
			cerr << "missing parameter for switch -" << TRpsbProcCmds::m_dimRpsbProcSwitches[TRpsbProcCmds::eOutFile].m_cToken << ", ignored.." << endl;
		else
			g_strDstFile = arg;
		break;
	case TRpsbProcCmds::eEVCutoff:
		if (NULL == arg)
			cerr << "missing parameter for switch -" << TRpsbProcCmds::m_dimRpsbProcSwitches[TRpsbProcCmds::eEVCutoff].m_cToken << ", ignored.." << endl;
		else
			g_dEValue = atof(arg);
		break;
	case TRpsbProcCmds::eQuiet:
		g_bSilent = true;
		break;
	case TRpsbProcCmds::eMode:
		if (NULL == arg)
			cerr << "missing parameter for switch -" << TRpsbProcCmds::m_dimRpsbProcSwitches[TRpsbProcCmds::eMode].m_cToken << ", ignored.." << endl;
		else
			g_iDataMode = GetIdx<TDataModes> (arg);
		break;
	case TRpsbProcCmds::eTData:	//target data
		if (NULL == arg)
			cerr << "missing parameter for switch -" << TRpsbProcCmds::m_dimRpsbProcSwitches[TRpsbProcCmds::eTData].m_cToken << ", ignored.." << endl;
		else
			g_iTargetData = GetIdx<TTargetData> (arg);
		break;
	case TRpsbProcCmds::eData:	//target data
		if (NULL == arg)
			cerr << "missing parameter for switch -" << TRpsbProcCmds::m_dimRpsbProcSwitches[TRpsbProcCmds::eData].m_cToken << ", ignored.." << endl;
		else
			g_strDataPath = arg;
		break;
	case TRpsbProcCmds::eFams:
		g_bSuperfams = true;
		break;
	case TRpsbProcCmds::eHelp:
		PrintUsage(g_pProgName, cout);
		g_bRunProg = false;
		break;
	case TRpsbProcCmds::eInvalid:
		cerr << "Invalid switch '" << arg[0] << "', ignored.." << endl;
		break;
	case TRpsbProcCmds::NONSWITCH:
	default:
		;
	}
}

/**********************************************************************
*	main
***********************************************************************/
int main(int argc, char * argv[])
{
	g_pProgName = argv[0];
	TRpsbProcCmds().IterParameters(argc, argv, ProcessCmdline);


	if (!g_bRunProg) return 2;
	
	if (!g_bSilent)	//display information
	{
		cerr << PROGRAM_TITLE << endl;
		cerr << "-------------------------------------------------------" << endl;
		cerr << "Download:" << endl;
		cerr << "\tftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/" << endl << endl;
		cerr << "Blast applications download:" << endl;
		cerr << "\tftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/" << endl << endl;
		cerr << "For more information please refer to the README file" << endl;
	}
	
	CNcbiRegistry reg;
	
	bool bHasConfig = false;
	
	if (!g_strCfgFile.empty())
		bHasConfig = ReadConfig(g_strCfgFile.c_str(), reg);
	else
		bHasConfig = ReadConfig((g_pProgName + string(".ini")).c_str(), reg) || ReadConfig(nullptr, reg);
	
	istream *pDataSrc = &cin;
	
	ifstream ifs;
	if (!g_strSrcFile.empty())
	{
		ifs.open(g_strSrcFile.c_str(), ios::binary | ios::in);
		pDataSrc = &ifs;
	}
	
	ostream *pDataDst = &cout;
	ofstream ofs;
	if (!g_strDstFile.empty())
	{
		ofs.open(g_strDstFile.c_str(), ios::binary | ios::out);
		pDataDst = &ofs;
	}
	// -- output to result
	(*pDataDst) << '#' << PROGRAM_TITLE << endl;
	
	(*pDataDst) << "#Config file:\t";
	if (bHasConfig)
	{
		if (g_strCfgFile.empty()) (*pDataDst) << g_pProgName << ".ini";
		else (*pDataDst) << g_strCfgFile;
	}
	(*pDataDst) << endl;
	
	(*pDataDst) << "#Input data file:\t";
	if (g_strSrcFile.empty()) (*pDataDst) << "stdin";
	else (*pDataDst) << g_strSrcFile;
	(*pDataDst) << endl;
	
	(*pDataDst) << "#Output data file:\t";
	if (g_strDstFile.empty()) (*pDataDst) << "stdout";
	else (*pDataDst) << g_strDstFile;
	(*pDataDst) << endl;
	(*pDataDst) << "#E-Value cutoff:\t" << g_dEValue << endl;
	(*pDataDst) << "#Redundancy:\t" << TDataModes::dimDisplay[g_iDataMode] << endl;
	(*pDataDst) << "#Data requested:\t" << TTargetData::dimDisplay[g_iTargetData] << endl;
	(*pDataDst) << "#Output format -- tab-delimited table" << endl;
	(*pDataDst) << '#' << DATASTART << endl;
	(*pDataDst) << '#' << SESSIONSTART << DELIMIT << "<session-ordinal>" << DELIMIT << "<program>" << DELIMIT << "<database>" << DELIMIT << "<score-matrix>" << DELIMIT << "<evalue-threshold>" << endl;
	(*pDataDst) << '#' << QUERYSTART << DELIMIT << "<query-id>" << DELIMIT << "<seq-type>" << DELIMIT << "<seq-length>" << DELIMIT << "<definition-line>" << endl;


	if (TTargetData::e_feats != g_iTargetData)
	{
		(*pDataDst) << '#' << DOMSTART << endl;
		(*pDataDst) << '#' << "<session-ordinal>" << DELIMIT << "<query-id[readingframe]>" << DELIMIT << "<hit-type>" << DELIMIT << "<PSSM-ID>" << DELIMIT << "<from>" << DELIMIT << "<to>" << DELIMIT << "<E-Value>" << DELIMIT << "<bitscore>" << DELIMIT << "<accession>" << DELIMIT << "<short-name>" << DELIMIT << "<incomplete>" << DELIMIT << "<superfamily PSSM-ID>" << endl;
		(*pDataDst) << "#more such lines......" << endl;
		(*pDataDst) << '#' << DOMEND << endl;
	}
	if (TTargetData::e_doms != g_iTargetData)
	{
		(*pDataDst) << '#' << FEATSTART << endl;
		(*pDataDst) << '#' << "<session-ordinal>" << DELIMIT << "<query-id[readingframe]>" << DELIMIT << "<annot-type>" << DELIMIT << "<title>" << DELIMIT << "<residue(coordinates)>" << DELIMIT << "<complete-size>" << DELIMIT << "<mapped-size>" << DELIMIT << "<source-domain>" << endl;
		(*pDataDst) << "#more such lines......" << endl;
		(*pDataDst) << '#' << FEATEND << endl;
		(*pDataDst) << '#' << MOTIFSTART << endl;
		(*pDataDst) << '#' << "<session-ordinal>" << DELIMIT << "<query-id[readingframe]>" << DELIMIT << "<title>" << DELIMIT << "<from>" << DELIMIT << "<to>" << DELIMIT << "<source-domain>" << endl;
		(*pDataDst) << "#more such lines......" << endl;
		(*pDataDst) << '#' << MOTIFEND << endl;
	}
	(*pDataDst) << '#' << QUERYEND << DELIMIT << "<query-id>" << endl;
	(*pDataDst) << "#more query sections.." << endl;
	(*pDataDst) << '#' << SESSIONEND << DELIMIT << "<session-ordinal>" << endl;
	(*pDataDst) << "#more session sections.." << endl;
	(*pDataDst) << '#' << DATAEND << endl;
	(*pDataDst) << "#=====================================================================" << endl;
	
	int objCount = 0;
	if (pDataSrc->good())
	{
		
		COflRpsbPostProcessor proc;
		if (bHasConfig)
			proc.LoadData(reg);
		else
			proc.LoadData();
		
		unsigned int procstatus = proc.GetStatus();
		
		if (procstatus > 0)
		{
			cerr << "Annotation data loading unsuccessful. Please check the follow data files to make sure they exist and are readable:" << endl;
			if (!g_strDataPath.empty())
			{
				if (procstatus & COflDomClstInfo::CDD_DATA_NOT_FOUND)
					cerr << '\t' << g_strDataPath << '/' << CDIDFILE << endl;
				if (procstatus & COflDomClstInfo::CLUSTER_LINK_NOT_FOUND)
					cerr << '\t' << g_strDataPath << '/' << CLSTLINKFILE << endl;
				if (procstatus & COflDomClstInfo::HIERARCHY_DATA_NOT_FOUND)
					cerr << '\t' << g_strDataPath << '/' << CDTRACKFILE << endl;
				if (procstatus & COflDomClstInfo::FEATURE_DATA_NOT_FOUND)
					cerr << '\t' << g_strDataPath << '/' << SPFEATFILE << endl;
				if (procstatus & COflDomClstInfo::GENERIC_FEATURE_DATA_NOT_FOUND)
					cerr << '\t' << g_strDataPath << '/' << GENFEATFILE << endl;
				if (procstatus & COflDomClstInfo::SPTHRESHOLD_DATA_NOT_FOUND)
					cerr << '\t' << g_strDataPath << '/' << MINBSCOREFILE << endl;
			}
			else if (bHasConfig)
			{
				if (procstatus & COflDomClstInfo::CDD_DATA_NOT_FOUND)
					cerr << '\t' << reg.Get(DATAPATH, CDDIDS) << endl;
				if (procstatus & COflDomClstInfo::CLUSTER_LINK_NOT_FOUND)
					cerr << '\t' << reg.Get(DATAPATH, CLUSTERLINKS) << endl;
				if (procstatus & COflDomClstInfo::HIERARCHY_DATA_NOT_FOUND)
					cerr << '\t' << reg.Get(DATAPATH, CDTRACKINFO) << endl;
				if (procstatus & COflDomClstInfo::FEATURE_DATA_NOT_FOUND)
					cerr << '\t' << reg.Get(DATAPATH, FEATURES) << endl;
				if (procstatus & COflDomClstInfo::GENERIC_FEATURE_DATA_NOT_FOUND)
					cerr << '\t' << reg.Get(DATAPATH, GENERIC_FEATURES) << endl;
				if (procstatus & COflDomClstInfo::SPTHRESHOLD_DATA_NOT_FOUND)
					cerr << '\t' << reg.Get(DATAPATH, SPECIFICTHRESHOLDS) << endl;
			}
			else
			{
				if (procstatus & COflDomClstInfo::CDD_DATA_NOT_FOUND)
					cerr << '\t' << CDIDFILE << endl;
				if (procstatus & COflDomClstInfo::CLUSTER_LINK_NOT_FOUND)
					cerr << '\t' << CLSTLINKFILE << endl;
				if (procstatus & COflDomClstInfo::HIERARCHY_DATA_NOT_FOUND)
					cerr << '\t' << CDTRACKFILE << endl;
				if (procstatus & COflDomClstInfo::FEATURE_DATA_NOT_FOUND)
					cerr << '\t' << SPFEATFILE << endl;
				if (procstatus & COflDomClstInfo::GENERIC_FEATURE_DATA_NOT_FOUND)
					cerr << '\t' << GENFEATFILE << endl;
				if (procstatus & COflDomClstInfo::SPTHRESHOLD_DATA_NOT_FOUND)
					cerr << '\t' << MINBSCOREFILE << endl;
			}
			return 255;
		}
		
		(*pDataDst) << DATASTART << endl;
		if (!g_bSilent) cerr << "Start process data..." << endl;
		objCount = proc.ProcessRpsbDataStream(*pDataSrc, *pDataDst, g_dEValue, g_iDataMode);
		if (!g_bSilent) cerr << "End process data..." << endl;
		(*pDataDst) << DATAEND << endl;
	}
	
	(*pDataDst) << endl << "#Total Blastout object processed\t" << objCount << endl;
	
	ifs.close();
	ofs.close();
	
	cerr << "Total CBlastOutput objects: " << objCount << endl;
	return 0;
}



