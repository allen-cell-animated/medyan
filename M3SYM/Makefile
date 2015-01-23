
#LIST ALL SOURCES HERE

SOURCES = $(wildcard *.cpp)
HEADERS = $(wildcard *.h)
OBJECTS = $(SOURCES:.cpp=.o)

#Configured for DEEPTHOUGHT2, UMD
CXX = /usr/bin/g++
DEBUG = -g
FAST = -O3 -ffast-math -mtune=corei7-avx -flto -funroll-loops -DNDEBUG -mtune=native

LDLIBS = -L/usr/local/lib/ -lboost_system -lgsl

#Make sure to add all folders for includes here
CPPFLAGS = -Wall -I/usr/local/include           \
		   -I./ -I./Structure -I./Chemistry     \
		   -I./Mechanics/Minimizer              \
		   -I./Mechanics/ForceField             \ 
		   -I./Mechanics/ForceField/Volume      \
		   -I./Mechanics/ForceField/Boundary    \
		   -I./Mechanics/ForceField/Branching   \
		   -I./Mechanics/ForceField/Filament    \
		   -I./Mechanics/ForceField/Linker      \
		   -I./Mechanics/ForceField/MotorGhost  \

LFLAGS = -Wall -std=c++11 -stdlib=libc++ $(FAST)

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
