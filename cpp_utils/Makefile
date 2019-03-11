CXX = g++

# compiler warning flags
WFLAGS = -Wall -W -Wextra -pedantic -Wno-attributes

# debug build flags
#DBFLAGS = -ggdb -pg
DBFLAGS = -ggdb -static
CXXFLAGS = $(WFLAGS) $(DBFLAGS)

# performance build flags
OFLAGS = -O3
CXXFLAGS = $(WFLAGS) $(OFLAGS)

# external libraries
LDLIBS = -lstdc++
#LDFLAGS = -pg

PROGS = table_parser\
	get_fasta_by_id	

PROGSOBJS = table_parser.o\
	get_fasta_by_id.o	

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
