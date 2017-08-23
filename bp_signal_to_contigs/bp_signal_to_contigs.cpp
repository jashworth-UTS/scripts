// parses samtools mpileup file, adds up counts on fwd and rvs strands, outputs these to separate files

//simple method to determine contiguous regions in position data (such as samtools depth)
//combines across gaps (-g)
//filters out small contigs (-m)

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <list>
#include <string>
#include <algorithm>
#include <cstdlib>

typedef std::string string;
typedef std::ifstream ifstream;
typedef std::ofstream ofstream;

void usage_error(){
	std::cerr << "\n"
	<< " -t #                 : threshold\n"
	<< " -g #                 : gap allowance\n"
	<< " -m #                 : min contig length\n"
	<< "\n";
	exit(EXIT_FAILURE);
}

typedef std::list<int> Positions;
typedef std::map<std::string, Positions> Seqs;
class Contig{
public:
	int start;
	int end;
	Contig(int s, int e){start=s; end=e;}
};
typedef std::list<Contig> Contigs;

int main(int argc, char* argv[]) {

	int threshold(0), gap(0), min(0);
	std::list<string> input;

	if(argc<2) usage_error();

	for (int i(1); i<argc; ++i) {
		string const & arg( argv[i] );

		if (arg == "-t") {
			if (++i >= argc) usage_error();
			threshold = atoi(argv[i]);

		} else if (arg == "-g") {
			if (++i >= argc) usage_error();
			gap = atoi(argv[i]);

		} else if (arg == "-m") {
			if (++i >= argc) usage_error();
			min = atoi(argv[i]);

		} else if (arg == "-h" || arg == "--help") usage_error();

		else input.push_back(arg);
	}

	Seqs seqs;

	if (input.size() < 1) usage_error();

	for (std::list<string>::const_iterator in(input.begin()); in != input.end(); ++in) {
		std::cerr << *in << std::endl;

		ifstream infile;
		infile.open(in->c_str());
		string line;
		while( getline(infile,line) ) {
			std::istringstream is(line);

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

	for(Seqs::iterator s(seqs.begin()); s!=seqs.end(); ++s){

		Contigs contigs;
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
				if (lastp-start >= min) contigs.push_back( Contig(start,lastp) );
				start = *p;
			}
			lastp = *p;
		}
		if (lastp-start >= min) contigs.push_back( Contig(start,lastp) );

		for(Contigs::iterator c(contigs.begin()); c!=contigs.end(); ++c){
//			std::cout << s->first << '\t' << c->start << '\t' << c->end << std::endl;
// pseudo-GTF
			std::cout << s->first <<'\t'<< "JA" <<'\t'<< "RNAseq_block" <<'\t'<< c->start <<'\t'<< c->end << "\t.\t+\t.\t" << "transcript_id \"" << s->first <<'_'<< c->start <<'_'<< c->end << "\"" << std::endl;
		}

	}

	return 0;
}
