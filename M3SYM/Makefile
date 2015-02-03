#LINUX MAKEFILE FOR M3SYM
SOURCES := $(shell find . -name '*.cpp')
HEADERS := $(shell find . -name '*.h')
OBJECTS := $(SOURCES:.cpp=.o)

#compiler choice and flags
CXX = /usr/bin/g++ -std=c++11
DEBUG = -g
FAST = -O3 -ffast-math -mtune=corei7-avx -flto -funroll-loops -DNDEBUG -mtune=native

LDLIBS = -L/usr/local/lib/ -lboost_system -lgsl 

CPPFLAGS = -Wall -Wno-sign-compare \
-I/usr/local/include -I/usr/include/boost \
-I./ -IStructure -IChemistry             \
-IMechanics/Minimizer                    \
-IMechanics/ForceField                   \
-IMechanics/ForceField/Volume            \
-IMechanics/ForceField/Boundary          \
-IMechanics/ForceField/Branching         \
-IMechanics/ForceField/Filament          \
-IMechanics/ForceField/Linker            \
-IMechanics/ForceField/MotorGhost        \

#MAIN USER MACROS
#Initialize chemical components of system
CPPFLAGS += -DCHEMISTRY
#Initialize mechanical components of system
CPPFLAGS += -DMECHANICS
#Use dynamic rate changing 
CPPFLAGS += -DDYNAMICRATES

#Boost memory macros
CPPFLAGS += -DBOOST_MEM_POOL
CPPFLAGS += -DBOOL_POOL_NSIZE=65536

#if compiling the testing suite
#CPPFLAGS += -DTESTING

#Detailed chemistry macros. 
#FOR ALL STANDARD FUNCTIONALITY, NONE OF THESE SHOULD BE TURNED OFF!!!

#Track dependents. Needed for Gillespie and NRM algorithm
CPPFLAGS += -DTRACK_DEPENDENTS

#Track zero copy and max copy number. 
#For passivating and activating reactions accordingly
CPPFLAGS += -DTRACK_ZERO_COPY_N
CPPFLAGS += -DTRACK_UPPER_COPY_N

#Reaction signaling. Needed for all Filament, Linker, 
#Motor, and BranchingPoint reactions in system
CPPFLAGS += -DREACTION_SIGNALING

#Species signaling
CPPFLAGS += -DRSPECIES_SIGNALING

LFLAGS = -Wall $(FAST)

all: M3SYM

M3SYM: $(OBJECTS)
	$(CXX) $(LDFLAGS) -o M3SYM $(OBJECTS) $(LDLIBS) $(FAST)

%.o: %.cpp 
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

clean:
	\find . -name '*.o' -delete
	\rm -f M3SYM

Makefile.dep: $(SOURCES) $(HEADERS)
	@echo "Updating the dependancy file: Makefile.dep"
	$(SHELL) -ec '$(CXX) -MM $(CPPFLAGS) $(INCLUDES)  $(filter %.cpp,$^) \
	| sed '\''s/^\(.*\)\.o[ :]*/\1.o $@ : /g'\'' > $@; \
	[ -s $@ ] || rm -f $@'

include Makefile.dep
