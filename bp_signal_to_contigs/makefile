CXX = g++

# compiler warning flags
WFLAGS = -Wall -W -Wextra -pedantic -std=c++11
DBFLAGS = -ggdb
OFLAGS = -O3

# external libraries
#LDLIBS = -lstdc++ -L/Users/jashwort/install_files/boost_1_56_0/stage/lib
#IFLAGS = -I/opt/local/include -I/Users/jashwort/install_files/boost_1_56_0

#CXXFLAGS = $(WFLAGS) $(DBFLAGS) $(IFLAGS)
CXXFLAGS = $(WFLAGS) $(OFLAGS) $(IFLAGS)

# build targets

all: bp_signal_to_contigs depth_to_utrs

bp_signal_to_contigs: bp_signal_to_contigs.cpp
	$(CXX) $(CXXFLAGS) -o bp_signal_to_contigs bp_signal_to_contigs.cpp

depth_to_utrs: depth_to_utrs.cpp
	$(CXX) $(CXXFLAGS) -o depth_to_utrs depth_to_utrs.cpp

clean:
	-rm -f *.o $(all)

.PHONY: tags
tags:
	ctags *.cpp *.h
