SOURCES = $(wildcard *.cpp)
HEADERS = $(wildcard *.h)
OBJECTS = $(SOURCES:.cpp=.o)

CXX = /opt/local/bin/g++-mp-4.7
DEBUG = -g
FAST = -O3
LDLIBS = -L/usr/local/lib/ -lboost_serialization 

# $(DEBUG) $(FAST)

CPPFLAGS = -Wall -I/usr/local/include -std=c++11 $(FAST)	 -c 
LFLAGS = -Wall -std=c++11 $(FAST)



all: CytoSim

CytoSim: $(OBJECTS)
	$(CXX) $(LDFLAGS) -o CytoSim $(OBJECTS) $(LDLIBS)

%.o: %.cpp 
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c $<

clean:
	\rm -f *.o
	\rm -f CytoSim

Makefile.dep: $(SOURCES) $(HEADERS)
	@echo "Updating the dependancy file: Makefile.dep"
	$(SHELL) -ec '$(CXX) -MM $(CPPFLAGS) $(INCLUDES)  $(filter %.cpp,$^) \
	| sed '\''s/^\(.*\)\.o[ :]*/\1.o $@ : /g'\'' > $@; \
	[ -s $@ ] || rm -f $@'


include Makefile.dep