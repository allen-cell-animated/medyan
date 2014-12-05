
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_common_h
#define M3SYM_common_h

#include <iostream>

#include "utility.h"

///Defining chemical and mechanical capablilites turned on/off
#define MECHANICS
//#define CHEMISTRY

///If compiling for testing
#define TESTING

///Other chemistry macros
#define TRACK_DEPENDENTS
#define TRACK_ZERO_COPY_N
#define TRACK_UPPER_COPY_N
#define RSPECIES_SIGNALING
#define REACTION_SIGNALING
//#define BOOST_MEM_POOL
//#define BOOL_POOL_NSIZE 65536

///Species constants
typedef unsigned short species_copy_t;
const species_copy_t max_ulim = 10000;
extern double global_time;

inline double tau() {return global_time;}

///To use STL containers, libraries, etc
using namespace std;

#endif
