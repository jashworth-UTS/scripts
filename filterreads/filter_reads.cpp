// JA 2018 this is high-speed/low-BS/brute force fastq read filtering script (by read IDs, e.g. hisat2 output)
// filterreads [fastqfile] [reads_to_exclude]
#include <map>
#include <string>
#include <list>
#include <iostream>
#include <fstream>
#include <cstdlib> // exit, EXIT_FAILURE

void help(){
	std::cerr << "usage: filterreads [fastqfile] [reads_to_exclude] (see example files)" << std::endl;
	std::cerr << "[reads_to_exclude] should be hisat2/sam output and/or have a ReadID\t... format" << std::endl;
}

typedef std::list<std::string> Strings;

int main(int argc, char *argv[]) {

	if (argc < 3){
		help();
		exit(EXIT_FAILURE);
	}

	std::string fastq_file, exclude_file;
	int i(1);
	fastq_file = argv[i++];
	exclude_file = argv[i++];

	std::map<std::string, bool> exclude;

	std::ifstream file;

	file.open( exclude_file.c_str() );
	if ( !file ) {
		std::cerr << "ERROR: unable to open file " << exclude_file.c_str() << std::endl;
		exit(EXIT_FAILURE);
	}
//	std::cerr << "Reading file " << exclude_file << std::endl;

	// read IDs to exclude (hisat2/sam format)
	std::string line;
	while ( getline(file,line) ){
		if(line[0]=='@') continue;
		std::string id(line.substr(0,line.find('\t')));
//		std::cerr << "ID line with ID: " << id << "." << std::endl;
//		std::cerr << id << std::endl;
		exclude[id] = true;
	}
	file.close();

	// "stream" filter fastq file
	file.open( fastq_file.c_str() );
	if ( !file ) {
		std::cerr << "ERROR: unable to open file " << fastq_file.c_str() << std::endl;
		exit(EXIT_FAILURE);
	}
//	std::cerr << "Reading file " << fastq_file << std::endl;

	// filter fastq reads
	// this flag is a simple toggler that switches on and off the printing of the line stream depending on whether the last ID encountered is in the exclude list or not
	bool filter(false);
	while ( getline(file,line) ){
		if(line[0]=='@'){
			std::string id(line.substr(1,line.find(' ')-1));
//			std::cerr << "ID line with ID: " << id << "." << std::endl;
			if (exclude.find(id) != exclude.end()) filter = true;
			else filter = false;
		}
		if(!filter) std::cout << line << '\n';
	}
	file.close();

}
