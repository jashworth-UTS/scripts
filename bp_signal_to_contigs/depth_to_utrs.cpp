// JA 2017
// read samtools depth
// create contigs above depth threshold (can span gaps, minimum lengths, see params)
// compare to genes in GTF file
// generate 'safe' putative 5 prime and 3 prime UTR extensions
// output in GTF format (for subsequent inclusion with input GTF)

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <list>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <regex>

typedef std::string string;
typedef std::ifstream ifstream;
typedef std::ofstream ofstream;
typedef std::istringstream istringstream;
typedef std::stringstream stringstream;

void usage_error(){
	std::cerr << "\n"
	<< " -t #                 : threshold\n"
	<< " -g #                 : gap allowance\n"
	<< " -m #                 : min contig length\n"
	<< " -gtf [file]          : gtf for filtering\n"
	<< " -r [regex]           : regex for gtf/gff geneids, e.g.:\n"
	<< "                      \"protein[idID_]*[= ]([A-Za-z0-9.-_]+)\"\n"
	<< "                      \"exon.*gene_id \\\"([A-Za-z0-9.-_]+)\\\"\"\n"
	<< "\n";
	exit(EXIT_FAILURE);
}

typedef std::list<int> Positions;
typedef std::list<string> Strings;
typedef std::map<std::string, Positions> SeqPositions;

class Contig{
public:
	string seq;
	bool strand; // is 'true' for sense/+/fwd strand
	int start;
	int end;
	string id;
	string parent;
	Contig() : seq(""), strand(true), start(0), end(0), id("") {}
	Contig(string _seq, bool _strand, int s, int e, string _id, string _parent="") : seq(_seq), strand(_strand), start(s), end(e), id(_id), parent(_parent) {}
	bool overlaps(Contig const &);
	bool contains (int pos) const { return (pos >= start && pos <= end); }
	bool contains (Contig const & other) const {
		// is other fully contained within self?
		return start<=other.start && end>=other.end;
	}
};

bool
Contig::overlaps(Contig const & other){
	// different scaffold (e.g. chromosome)
	if(seq!=other.seq) return false;
	// is fully left of other
	if(start<other.start && end<other.start) return false;
	// is fully right of other
	if(start>other.end && end>other.end) return false;
	// is fully contained within other
	if(start>other.start && end<other.end) return false;
	return true;
}

bool operator < (Contig const & c1, Contig const & c2){
if(c1.seq<c2.seq) return true;
if(c1.seq>c2.seq) return false;
if(c1.start<c2.start) return true;
if(c1.start>c2.start) return false;
if(c1.end<c2.end) return true;
if(c1.end>c2.end) return false;
return false;
}

typedef std::list<Contig> Contigs;
typedef std::map<string, Contigs> SeqContigs;

// this extra derived class is a bit superfluous at the moment
class Gene : public Contig {
public:
	void dimension(string _seq, int _start, int _end, bool _strand){
		if(_seq != seq) { std::cerr << "WARNING: seq mismatch for gene id " << id << std::endl; return; }
		if(_strand != strand) { std::cerr << "WARNING: strand mismatch for gene id " << id << std::endl; return; }
		if(_start < start) start = _start;
		if(_end > end) end = _end;
	}
	Gene() : Contig() {}
	Gene(string seq, bool strand, int start, int end, string id) :
		Contig(seq,strand,start,end,id){}
};

typedef std::map<string,Gene> Genes;
typedef std::map<string,Genes> Genome;

// this could be made into a << operator for SeqContigs
void print_contigs(SeqContigs const & contigs, string type, bool gtf=true){
	for(SeqContigs::const_iterator sc(contigs.begin()); sc!=contigs.end(); ++sc){
		Contigs const & ctgs(contigs.at(sc->first));
		for(Contigs::const_iterator c(ctgs.begin()); c!=ctgs.end(); ++c){
			string strand("+");
			if(!c->strand) strand="-";
			if(gtf) std::cout << c->seq <<'\t'<< "JA" <<'\t'<< type << '\t'<< c->start <<'\t'<< c->end << "\t.\t" << strand << "\t.\t" << "parent \"" << c->parent << "\"; "<< "id \"" << c->id << "\"; " << std::endl;
			else std::cout << type << '\t' << c->seq << '\t' << c->start << '\t' << c->end << '\t' << c->id << std::endl;
		}
	}
}

int main(int argc, char* argv[]) {

	int threshold(0), gap(0), min(0);
	string gtffile, geneid_regex;
	std::list<string> input;

	// for Thaps3All or Thaps3Filtered gffs
	geneid_regex = "protein[idID_]*[= ]([A-Za-z0-9.-_]+)";
	// for Braker1/Augustus as of Aug '17
//	geneid_regex = "exon.*gene_id \"([A-Za-z0-9.-_]+)\"";

	if(argc<2) usage_error();

	for (int i(1); i<argc; ++i) {
		string const & arg( argv[i] );

		if (arg == "-t") {
			if (++i >= argc) usage_error();
			threshold = atoi(argv[i]);

		} else if (arg == "-g") {
			if (++i >= argc) usage_error();
			gap = atoi(argv[i]);

		} else if (arg == "-r") {
			if (++i >= argc) usage_error();
			geneid_regex = argv[i];

		} else if (arg == "-gtf") {
			if (++i >= argc) usage_error();
			gtffile = argv[i];

		} else if (arg == "-m") {
			if (++i >= argc) usage_error();
			min = atoi(argv[i]);

		} else if (arg == "-h" || arg == "--help") usage_error();

		else input.push_back(arg);
	}

	SeqPositions seqs;

	if (input.size() < 1) usage_error();

	for (std::list<string>::const_iterator in(input.begin()); in != input.end(); ++in) {
		std::cerr << *in << std::endl;

		ifstream infile;
		infile.open(in->c_str());
		string line;
		while( getline(infile,line) ) {
			istringstream is(line);

			// get first three fields
			string seqn; int pos; int depth;
			is >> seqn;
			is >> pos;
			is >> depth;

			if(depth < threshold) continue;

			if(seqs.find(seqn)==seqs.end()) seqs[seqn] = Positions();

			seqs[seqn].push_back(pos);
		}
	}

	SeqContigs contigs;

	for(SeqPositions::iterator s(seqs.begin()); s!=seqs.end(); ++s){

		contigs[s->first] = Contigs();
		Contigs & ctgs(contigs[s->first]);
		bool first(true);
		int start(0);
		int lastp(0);

		Positions & ps(s->second);
		ps.sort();
		for(Positions::iterator p(ps.begin()); p!=ps.end(); ++p){
			if(first){
				start = *p;
				first = false;
			} else if (*p > lastp + gap + 1){
				if (lastp-start >= min) ctgs.push_back( Contig(s->first,true,start,lastp,"") );
				start = *p;
			}
			lastp = *p;
		}
		if (lastp-start >= min) ctgs.push_back( Contig(s->first,true,start,lastp,"") );
	}

//	stringstream ss;
//	ss << "RNAseq_t" << threshold << "_m" << min;
//	print_contigs(contigs,ss.str());

	// filter down to contigs that could be UTRs
	// filter by overlap with first exon/CDS start, or last exon/CDS end

	// load annotated genes from GTF/GFF

	Genome genome;

	std::regex gidrgx(geneid_regex);

	ifstream infile;
	infile.open(gtffile.c_str());
	string line;
	while( getline(infile,line) ) {
		std::smatch match;
		if(!std::regex_search(line,match,gidrgx)) continue;
		istringstream is(line);
		string seqn, source, type, score, strand, frame;
		int start, end;
		is >> seqn >> source >> type >> start >> end >> score >> strand >> frame;
		if(genome.find(seqn)==genome.end()) genome[seqn] = Genes();
		Genes & genes(genome[seqn]);
		string id(match.str(1));
		bool strandbool(true);
		if(strand!="+") strandbool=false;
		if(genes.find(id)==genes.end()) genes[id] = Gene(seqn,strandbool,start,end,id);
		else genes[id].dimension(seqn,start,end,strandbool);
	}

	SeqContigs utrs_upstream;
	SeqContigs utrs_downstream;

	// filter for contigs overlapping with ends of only one gene, trim to gene boudaries, rename to reflect putative parent and UTR position
	for(SeqContigs::iterator sc(contigs.begin()); sc!=contigs.end(); ++sc){
		string seqn(sc->first);
		utrs_upstream[seqn] = Contigs();
		utrs_downstream[seqn] = Contigs();
		Genes & genes(genome[seqn]);
		Contigs & ctgs(contigs[seqn]);
		stringstream ss;
		for(Contigs::iterator c(ctgs.begin()); c!=ctgs.end(); ++c){
			typedef std::list<Genes::const_iterator> GeneIts;
			GeneIts parents;
			// overlap?
			for(Genes::iterator g(genes.begin()); g!=genes.end(); ++g){
				if(c->overlaps(g->second)) parents.push_back(g);
			}

			// single overlap?
			// unfortunately too simple: lack of directionality logic results in false negatives
			//if(parents.size()!=1) continue;

			// strand logic, trim, rename and store
			for(GeneIts::iterator g1(parents.begin()); g1!=parents.end(); ++g1){
				Gene const & parent((*g1)->second);
				bool strand(parent.strand);
				// start is free from g1 context
				if(!parent.contains(c->start)){
					// check other parents for valid start
					bool valid(true);
					for(GeneIts::iterator g2(parents.begin()); g2!=parents.end(); ++g2){
						Gene const & other((*g2)->second);
						if(g1==g2) continue;
						if(other.contains(c->start)) {valid=false; break;}
						if(other.start<parent.start) {valid=false; break;}
					}
					if(valid){
						if(strand){
							ss.str("");
							ss << parent.id << "_" << "5pUTR_RNAseq";
							utrs_upstream[seqn].push_back( Contig(seqn,strand,c->start,parent.start-1,ss.str(),parent.id) );
						} else {
							ss.str("");
							ss << parent.id << "_" << "3pUTR_RNAseq";
							utrs_downstream[seqn].push_back( Contig(seqn,strand,c->start,parent.start-1,ss.str(),parent.id) );
						}
					}
				}

				if(!parent.contains(c->end)){
					// check other parents for valid end
					bool valid(true);
					for(GeneIts::iterator g2(parents.begin()); g2!=parents.end(); ++g2){
						Gene const & other((*g2)->second);
						if(g1==g2) continue;
						if(other.contains(c->end)) {valid=false; break;}
						if(other.end>parent.end) {valid=false; break;}
					}
					if(valid){
						if(strand) {
							ss.str("");
							ss << parent.id << "_" << "3pUTR_RNAseq";
							utrs_downstream[seqn].push_back( Contig(seqn,strand,parent.end+1,c->end,ss.str(),parent.id) );
						} else {
							ss.str("");
							ss << parent.id << "_" << "5pUTR_RNAseq";
							utrs_upstream[seqn].push_back( Contig(seqn,strand,parent.end+1,c->end,ss.str(),parent.id) );
						}
					}
				}
			}
		}
	}

	stringstream ss;
	ss << "five_prime_utr";
	print_contigs(utrs_upstream,ss.str());

	ss.str("");
	ss << "three_prime_utr";
	print_contigs(utrs_downstream,ss.str());

	return 0;
}
