SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

CXX = /opt/local/bin/g++-mp-4.7
DEBUG = -g
FAST = -O3 -ffast-math
LDLIBS = -L../ -L/usr/local/lib/ -L/usr/local/gcc47/gtest-1.6.0/lib/.libs -lcyto -lgtest -lgtest_main -lboost_system

# $(DEBUG) $(FAST)

CPPFLAGS = -Wall -I../ -I/usr/local/include -I/usr/local/gcc47/gtest-1.6.0/include -std=c++11 $(FAST) 
LFLAGS = -Wall -std=c++11 $(FAST)

all: gtest-cyto test

gtest-cyto: $(OBJECTS)
	$(CXX) $(LDFLAGS) -o gtest-cyto $(OBJECTS) $(LDLIBS)

test: gtest-cyto
	export DYLD_LIBRARY_PATH=/usr/local/gcc47/gtest-1.6.0/lib/.libs:/Users/papoian/Dropbox/Codes/CytoSim/CytoSim:/usr/local/lib/:/opt/local/lib
	./gtest-cyto

%.o: %.cpp 
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c $<

clean:
	\rm -f *.o
	\rm -f gtest-cyto

Makefile.dep: $(SOURCES) $(HEADERS)
	@echo "Updating the dependancy file: Makefile.dep"
	$(SHELL) -ec '$(CXX) -MM $(CPPFLAGS) $(INCLUDES)  $(filter %.cpp,$^) \
	| sed '\''s/^\(.*\)\.o[ :]*/\1.o $@ : /g'\'' > $@; \
	[ -s $@ ] || rm -f $@'

include Makefile.dep