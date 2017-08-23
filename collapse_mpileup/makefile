CXX = g++

# compiler warning flags
WFLAGS = -Wall -W -Wextra -pedantic
DBFLAGS = -ggdb
OFLAGS = -O3
IFLAGS = -I/opt/local/include -I/Users/jashwort/install_files/boost_1_56_0
#CXXFLAGS = $(WFLAGS) $(DBFLAGS) $(IFLAGS)
CXXFLAGS = $(WFLAGS) $(OFLAGS) $(IFLAGS)

EXE = collapse_mpileup.bin
OBJECTFILES = collapse_mpileup.o

# external libraries
LDLIBS = -lstdc++ -L/Users/jashwort/install_files/boost_1_56_0/stage/lib

# build targets
all: $(EXE)
$(EXE): $(OBJECTFILES)
	$(CXX) $(OBJECTFILES) $(LDLIBS) -o $(EXE)

clean:
	-rm *.o $(EXE)

.PHONY: tags
tags:
	ctags *.cpp *.h
