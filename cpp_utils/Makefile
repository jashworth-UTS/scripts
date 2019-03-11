CXX = g++

# compiler warning flags
WFLAGS = -Wall -W -Wextra -pedantic -Wno-attributes

# debug build flags
#DBFLAGS = -ggdb -pg
DBFLAGS = -ggdb -static
CMPLFLAGS = -std=c++11
CXXFLAGS = $(WFLAGS) $(DBFLAGS) $(CMPLFLAGS)

# performance build flags
OFLAGS = -O3
CXXFLAGS = $(WFLAGS) $(OFLAGS) $(CMPLFLAGS)

# external libraries
LDLIBS = -lstdc++
#LDFLAGS = -pg

PROGS = table_parser\
	get_fasta_by_id\
	filter_reads

PROGSOBJS = table_parser.o\
	get_fasta_by_id.o\
	filter_reads.o

#EXE = get_fasta_by_id
#OBJECTFILES = get_fasta_by_id.o

# build targets
#all: $(EXE)
all: $(PROGS)

clean:
	-rm ${PROGS} ${PROGOBJS}

#$(EXE): $(OBJECTFILES)
##	$(CXX) $(OBJECTFILES) $(LDLIBS) -o $(EXE) $(LDFLAGS)
#	$(CXX) $(OBJECTFILES) $(LDLIBS) -o $(EXE)

.PHONY: tags
tags:
	ctags *.cpp *.h
