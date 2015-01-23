#LINUX MAKEFILE FOR M3SYM
#Configured for DEEPTHOUGHT2, UMD

#LIST ALL SOURCE DIRECTORIES HERE
VPATH = Structure:Chemistry:Mechanics/Minimizer:Mechanics/ForceField: \
Mechanics/ForceField/Volume:Mechanics/ForceField/Boundary:            \
Mechanics/ForceField/Branching:Mechanics/ForceField/Filament:	      \
Mechanics/ForceField/Linker:Mechanics/ForceField/MotorGhost

SOURCES = $(wildcard *.cpp)
HEADERS = $(wildcard *.h)
OBJECTS = $(SOURCES:.cpp=.o)

CXX = /usr/bin/c++ -std=c++0x
DEBUG = -g
FAST = -O3 -ffast-math -mtune=corei7-avx -flto -funroll-loops -DNDEBUG -mtune=native

LDLIBS = -L/usr/local/lib/ -lboost_system -lgsl 

#Make sure to add all folders for includes here
CPPFLAGS = -Wall -I/usr/local/include    \
-I/usr/include/boost                     \
-I./ -IStructure -IChemistry             \
-IMechanics/Minimizer                    \
-IMechanics/ForceField                   \
-IMechanics/ForceField/Volume            \
-IMechanics/ForceField/Boundary          \
-IMechanics/ForceField/Branching         \
-IMechanics/ForceField/Filament          \
-IMechanics/ForceField/Linker            \
-IMechanics/ForceField/MotorGhost        \

LFLAGS = -Wall $(FAST)

all: M3SYM

M3SYM: $(OBJECTS)
	$(CXX) $(LDFLAGS) -o M3SYM $(OBJECTS) $(LDLIBS) $(FAST)

%.o: %.cpp 
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c $<

clean:
	\rm -f *.o
	\rm -f M3SYM

Makefile.dep: $(SOURCES) $(HEADERS)
	@echo "Updating the dependancy file: Makefile.dep"
	$(SHELL) -ec '$(CXX) -MM $(CPPFLAGS) $(INCLUDES)  $(filter %.cpp,$^) \
	| sed '\''s/^\(.*\)\.o[ :]*/\1.o $@ : /g'\'' > $@; \
	[ -s $@ ] || rm -f $@'

#include Makefile.dep
