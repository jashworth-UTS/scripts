// script for rapidly parsing large tsv tables (currently blastx outfmt 6 tables)
// c++ is faster than anything else for ripping through very large files

#include <unordered_map>
#include <set>
#include <string>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib> // exit, EXIT_FAILURE
#include <regex>
#include <thread>

void help(){
	std::cerr << "usage: table_parser tablefiles" << std::endl;
	exit(EXIT_FAILURE);
}

typedef std::thread thread;
typedef std::vector<thread> Threads;
typedef std::regex regex;
typedef std::string string;
typedef std::vector<string> Strings;
typedef std::set<string> StringSet;
typedef std::unordered_map<string, unsigned> Counts;
typedef std::unordered_map<string, Counts> Stat;
typedef std::unordered_map<string, Stat> Stats;

void output_json(Stats const & stats, std::ostream & out)
{
	out << "{";
	for(Stats::const_iterator fl(stats.begin()); fl!=stats.end(); ++fl){
		if(fl != stats.begin()) out << ",";
		out << "\"" << fl->first << "\"" << ":{";
		Stat const & flstats(fl->second);
		for(Stat::const_iterator fs(flstats.begin()); fs!=flstats.end(); ++fs){
			if(fs != flstats.begin()) out << ",";
			out << "\"" << fs->first << "\"" << ":{";
			Counts const & cnts(fs->second);
			for(Counts::const_iterator cnt(cnts.begin()); cnt!=cnts.end(); ++cnt){
				string key(cnt->first);
				key = std::regex_replace(key, std::regex("\\\\"), "");
				if(cnt != cnts.begin()) out << ",";
				out << "\"" << key << "\"" << ":" << cnt->second;
			}
			out << "}";
//			out << std::endl;
		}
		out << "}";
//		out << std::endl;

	}
	out << "}" << std::endl;
}

void output_tables(Stats const & stats, std::string sep="\t"){
	// collect all unique stat types
	StringSet types;
	// ordered set of sample names (important for ensured header/rows correspondence)
	StringSet samples;
	for(Stats::const_iterator fl(stats.begin()); fl!=stats.end(); ++fl){
		samples.insert(fl->first);
		Stat const & flstats(fl->second);
		for(Stat::const_iterator fs(flstats.begin()); fs!=flstats.end(); ++fs){
			types.insert(fs->first);
		}
	}

	std::copy(types.begin(), types.end(), std::ostream_iterator<std::string>(std::cerr, " "));

	// separate table files for each stat type
	for(StringSet::const_iterator tp(types.begin()); tp!=types.end(); ++tp){

		// collect all unique row ids for type
		StringSet allids;
		for(StringSet::const_iterator smp(samples.begin()); smp!=samples.end(); ++smp){
			Stat const & stat(stats.at(*smp));
			Counts const & cnts(stat.at(*tp));
			for(Counts::const_iterator cnt(cnts.begin()); cnt!=cnts.end(); ++cnt){
				allids.insert(cnt->first);
			}
		}

		std::ofstream ofs;
		std::string fname( (*tp) + ".tsv");
		ofs.open(fname.c_str());

		// header
		ofs << "id";
		for(StringSet::const_iterator smp(samples.begin()); smp!=samples.end(); ++smp){
			ofs << sep;
			ofs << *smp;
		}
		ofs << "\n";

		// output full row for each unique id
		for(StringSet::const_iterator id(allids.begin()); id!=allids.end(); ++id){
			ofs << *id;
			for(StringSet::const_iterator smp(samples.begin()); smp!=samples.end(); ++smp){
				ofs << sep;
				Stat const & stat(stats.at(*smp));
				Counts const & cnts(stat.at(*tp));
				if(cnts.find(*id) != cnts.end()) ofs << cnts.at(*id);
				else ofs << 0;
			}
			ofs << "\n";
		}

		ofs.close();
	}

}

void parse_table(string const & fname, Stat & filestats, regex const & rgx){

		std::ifstream file_o;
		file_o.open( fname.c_str() );
		if ( !file_o ) {
			std::cerr << "ERROR: unable to open file " << fname << std::endl;
			exit(EXIT_FAILURE);
		}
//		std::cerr << "Reading file " << fname << std::endl;

//		Counts incl;
		string line;
		while ( getline(file_o,line) ){
		// for tables with only one (top) hit (diamond: --max-target-seqs 1), skip id accounting, as it results in much larger memory footprint
//			string id(line.substr(0,line.find('\t')));
//			std::cerr << id;
//			if(incl.find(id) == incl.end()) incl[id] = 0;
//			else continue;
//			incl[id]++;
//			if(incl[id] > max_incl) continue;
			std::smatch matches;
			if(std::regex_search(line, matches, rgx)){
//				string acc(matches[0].str());
				string acc(matches[1].str());
				string os(matches[2].str());
				string ox(matches[3].str());
				string gn(matches[4].str());

//				// due to computational expense of re-parsing names millions of times, do this name collapsing later, post-tabulation (e.g. in Python scripts)
//				if(collapse_strains){
//					std::istringstream is(os);
//					os.clear();
//					string word;
//					unsigned n(0);
//					while(getline(is,word,' ') && n<2){ ++n; os += " " + word; }
//				}

//				std::cerr << acc << std::endl;
//				std::cerr << " " << os << " " << ox << " " << gn;
				if(filestats["acc"].find(acc) == filestats["acc"].end()) filestats["acc"][acc] = 0;
				if(filestats["os"].find(os) == filestats["os"].end()) filestats["os"][os] = 0;
				if(filestats["ox"].find(ox) == filestats["ox"].end()) filestats["ox"][ox] = 0;
				if(filestats["gn"].find(gn) == filestats["gn"].end()) filestats["gn"][gn] = 0;
				filestats["acc"][acc]++;
				filestats["os"][os]++;
				filestats["ox"][ox]++;
				filestats["gn"][gn]++;
			}
//			std::cerr << std::endl;
		}
		file_o.close();
	}

int main(int argc, char *argv[]) {

	if (argc < 2) help();

	Strings table_files;
	for(int i(1); i<argc; ++i){
		table_files.push_back(argv[i]);
	}

	size_t nthreads(table_files.size());

	std::cerr << "parsing files with " << nthreads << " threads:" << std::endl;
	std::copy(table_files.begin(), table_files.end(), std::ostream_iterator<string>(std::cerr, "\n"));
	std::cerr << std::flush;

//	unsigned max_incl(1);
//	bool collapse_strains(true);

	Stats stats;
	Strings stat_types;
	stat_types.push_back("acc");
	stat_types.push_back("os");
	stat_types.push_back("ox");
	stat_types.push_back("gn");

// A00152:58:H5KW2DMXX:2:1101:31241:1016	tr|K4FSH1|K4FSH1_CALMI	1.1e-08	66.6	33	97.0	100	2	tr|K4FSH1|K4FSH1_CALMI Tubulin beta chain OS=Callorhinchus milii OX=7868 PE=2 SV=1
	regex rgx("tr\\|([^\\|]+?)\\|.*?OS=(.+?) OX=(.+?) GN=([^ ]+)");

	Threads threads;
	for(unsigned i(0); i<table_files.size(); ++i){
		std::string const & fname(table_files[i]);
		// set up containers
		stats[fname] = Stat();
		Stat & filestats(stats[fname]);
		for(Strings::const_iterator type(stat_types.begin()); type!=stat_types.end(); ++type) filestats[*type] = Counts();
		threads.push_back( thread(parse_table, fname, std::ref(filestats), rgx) );
//		parse_table(fname, filestats, rgx);
	}

	std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));
	std::cerr << "threads done" << std::endl;

	output_json(stats,std::cout);
	output_tables(stats);
}
