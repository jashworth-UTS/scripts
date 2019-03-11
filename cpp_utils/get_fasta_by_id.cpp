// JA 2019 this is high-speed/low-BS/brute force fasta filtering script (by fasta IDs)
// filterreads [ids] [fastafile]
#include <map>
#include <string>
#include <list>
#include <iostream>
#include <fstream>
#include <cstdlib> // exit, EXIT_FAILURE

void help(){
	std::cerr << "usage: get_fasta_by_id [ids] [fastafile]" << std::endl;
}

typedef std::list<std::string> Strings;

int main(int argc, char *argv[]) {

	if (argc < 3){
		help();
		exit(EXIT_FAILURE);
	}

	std::string ids_file, fasta_file;
	int i(1);
	ids_file = argv[i++];
	fasta_file = argv[i++];

	std::map<std::string, bool> ids;

	std::ifstream file;

	file.open( ids_file.c_str() );
	if ( !file ) {
		std::cerr << "ERROR: unable to open file " << ids_file.c_str() << std::endl;
		exit(EXIT_FAILURE);
	}
//	std::cerr << "Reading file " << ids_file << std::endl;

	// read IDs to include
	std::string line;
	while ( getline(file,line) ){
		std::string id(line);
//		std::cerr << "ID line with ID: " << id << "." << std::endl;
//		std::cerr << id << std::endl;
		ids[id] = true;
	}
	file.close();

	// "stream" filter fasta file
	file.open( fasta_file.c_str() );
	if ( !file ) {
		std::cerr << "ERROR: unable to open file " << fasta_file.c_str() << std::endl;
		exit(EXIT_FAILURE);
	}
//	std::cerr << "Reading file " << fasta_file << std::endl;

	bool output(false);
	std::string id("");
	while ( getline(file,line) ){
		if(line[0]=='>'){
			id = line.substr(1,line.find(' ')-1);
//			std::cerr << "ID line with ID: " << id << "." << std::endl;
			if (ids.find(id) != ids.end()) output = true;
			else output = false;
		}
		if(output) std::cout << line << '\n';
	}
	file.close();

}
