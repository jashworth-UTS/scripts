#include <unordered_map>
#include <string>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib> // exit, EXIT_FAILURE
#include <regex>

void help(){
	std::cerr << "usage: table_parser [tablefile]" << std::endl;
}

typedef std::string string;
typedef std::list<string> Strings;
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

int main(int argc, char *argv[]) {

	if (argc < 2){
		help();
		exit(EXIT_FAILURE);
	}

	Strings table_files;
	for(int i(1); i<argc; ++i){
		table_files.push_back(argv[i]);
	}

	unsigned max_incl(1);
	Stats stats;
	Strings stat_types;
	stat_types.push_back("os");
	stat_types.push_back("ox");
	stat_types.push_back("gn");

//	std::regex rgx("tr\|([^\|]+?)\|.*?OS=(.*?) OX=(.*?) GN=([^ ]+)");
	std::regex rgx("OS=(.+?) OX=(.+?) GN=([^ ]+)");

	for(Strings::const_iterator file(table_files.begin()); file!=table_files.end(); ++file){
		stats[*file] = Stat();
		Stat & filestats(stats[*file]);
		for(Strings::const_iterator type(stat_types.begin()); type!=stat_types.end(); ++type){
			stats[*file][*type] = Counts();
		}
		Counts incl;

		std::ifstream file_o;
		file_o.open( file->c_str() );
		if ( !file_o ) {
			std::cerr << "ERROR: unable to open file " << file->c_str() << std::endl;
			exit(EXIT_FAILURE);
		}
//		std::cerr << "Reading file " << *file << std::endl;

		string line;
		while ( getline(file_o,line) ){
			std::istringstream linestream(line);
			string id;
			linestream >> id;
//			std::cerr << id;
			if(incl.find(id) == incl.end()) incl[id] = 0;
			incl[id]++;
			if(incl[id] > max_incl) continue;
			std::smatch matches;
			if(std::regex_search(line, matches, rgx)){
//				string acc(matches[0].str());
				string os(matches[1].str());
				string ox(matches[2].str());
				string gn(matches[3].str());
//				std::cerr << " " << os << " " << ox << " " << gn;
				if(filestats["os"].find(os) == filestats["os"].end()) filestats["os"][os] = 0;
				if(filestats["ox"].find(ox) == filestats["ox"].end()) filestats["ox"][ox] = 0;
				if(filestats["gn"].find(gn) == filestats["gn"].end()) filestats["gn"][gn] = 0;
				filestats["os"][os]++;
				filestats["ox"][ox]++;
				filestats["gn"][gn]++;
			}
//			std::cerr << std::endl;
		}
		file_o.close();
	}
	
	output_json(stats,std::cout);
}
