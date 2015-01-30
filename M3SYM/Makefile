#LINUX MAKEFILE FOR M3SYM
SOURCES := $(shell find . -name '*.cpp')
HEADERS := $(shell find . -name '*.h')
OBJECTS := $(SOURCES:.cpp=.o)

CXX = /usr/bin/g++ -std=c++11
DEBUG = -g
FAST = -O3 -ffast-math -mtune=corei7-avx -flto -funroll-loops -DNDEBUG -mtune=native

LDLIBS = -L/usr/local/lib/ -lboost_system -lgsl

CPPFLAGS = -Wall -I/usr/local/include -I/usr/include/boost \
-I./ -IStructure -IChemistry             \
-IMechanics/Minimizer                    \
-IMechanics/ForceField                   \
-IMechanics/ForceField/Volume            \
-IMechanics/ForceField/Boundary          \
-IMechanics/ForceField/Branching         \
-IMechanics/ForceField/Filament          \
-IMechanics/ForceField/Linker            \
-IMechanics/ForceField/MotorGhost        \

LFLAGS = $(FAST)

all: M3SYM

M3SYM: $(OBJECTS)
	$(CXX) $(LDFLAGS) -o M3SYM $(OBJECTS) $(LDLIBS) $(FAST)

%.o: %.cpp 
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

clean:
	\find . -name '*.o' -delete
	\rm -f M3SYM

#Makefile.dep: $(SOURCES) $(HEADERS)
#	@echo "Updating the dependancy file: Makefile.dep"
#	$(SHELL) -ec '$(CXX) -MM $(CPPFLAGS) $(INCLUDES)  $(filter %.cpp,$^) \
#	| sed '\''s/^\(.*\)\.o[ :]*/\1.o $@ : /g'\'' > $@; \
#	[ -s $@ ] || rm -f $@'

#include Makefile.dep
