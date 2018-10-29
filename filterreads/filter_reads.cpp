#include <map>
#include <string>
#include <list>
#include <iostream>
#include <fstream>
#include <cstdlib> // exit, EXIT_FAILURE

typedef std::list<std::string> Strings;

int main(int argc, char *argv[]) {

	std::string fqf, exf;
	int i(1);
	fqf = argv[i++];
	exf = argv[i++];

	std::map<std::string, bool> ex;

	std::ifstream file;

	file.open( exf.c_str() );
	if ( !file ) {
		std::cerr << "ERROR: unable to open file " << exf.c_str() << std::endl;
		exit(EXIT_FAILURE);
	}
//	std::cerr << "Reading file " << exf << std::endl;

	// read ids to exclude (sam format)
	std::string line;
	while ( getline(file,line) ){
		if(line[0]=='@') continue;
		std::string id(line.substr(0,line.find('\t')));
//		std::cerr << "ID line with ID: " << id << "." << std::endl;
//		std::cerr << id << std::endl;
		ex[id] = true;
	}
	file.close();

	// "stream" filter fastq file
	file.open( fqf.c_str() );
	if ( !file ) {
		std::cerr << "ERROR: unable to open file " << fqf.c_str() << std::endl;
		exit(EXIT_FAILURE);
	}
//	std::cerr << "Reading file " << fqf << std::endl;

	// read ids to exclude (sam format)
	bool filter(false);
	while ( getline(file,line) ){
		if(line[0]=='@'){
			std::string id(line.substr(1,line.find(' ')-1));
//			std::cerr << "ID line with ID: " << id << "." << std::endl;
			if (ex.find(id) != ex.end()) filter = true;
			else filter = false;
		}
		if(!filter) std::cout << line << '\n';
	}
	file.close();

}
