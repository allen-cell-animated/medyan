SOURCES = $(wildcard *.cpp)
HEADERS = $(wildcard *.h)
OBJECTS = $(SOURCES:.cpp=.o)

CXX = /opt/local/bin/g++-mp-4.7
DEBUG = -g
FAST = -O3 -ffast-math
LDLIBS = -L/usr/local/lib/ -lboost_serialization 

# $(DEBUG) $(FAST)

CPPFLAGS = -Wall -I/usr/local/include -std=c++11 $(FAST)
LFLAGS = -Wall -std=c++11 $(FAST)



all: CytoSim

CytoSim: $(OBJECTS)
	$(CXX) $(LDFLAGS) -o CytoSim $(OBJECTS) $(LDLIBS)

%.o: %.cpp 
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c $<

lib: ChemGillespieImpl.o ChemNRMImpl.o  ChemSim.o ChemSimpleGillespieImpl.o Component.o Compartment.o Composite.o  RSpecies.o Reaction.o Species.o System.o common.o
	$(CXX) $(CFLAGS) $(CPPFLAGS) -dynamiclib -undefined suppress -flat_namespace ChemGillespieImpl.o ChemNRMImpl.o  ChemSim.o ChemSimpleGillespieImpl.o Component.o Compartment.o Composite.o  RSpecies.o Reaction.o Species.o System.o common.o -o libcyto.dylib

clean:
	\rm -f *.o
	\rm -f CytoSim

Makefile.dep: $(SOURCES) $(HEADERS)
	@echo "Updating the dependancy file: Makefile.dep"
	$(SHELL) -ec '$(CXX) -MM $(CPPFLAGS) $(INCLUDES)  $(filter %.cpp,$^) \
	| sed '\''s/^\(.*\)\.o[ :]*/\1.o $@ : /g'\'' > $@; \
	[ -s $@ ] || rm -f $@'


include Makefile.dep