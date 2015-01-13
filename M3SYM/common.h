
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
#define CHEMISTRY

///Defining dynamic reaction rates (chemistry and mechanics both on)
#define DYNAMICRATES

///If compiling for testing
#define TESTING

///Other chemistry macros
#define TRACK_DEPENDENTS
#define TRACK_ZERO_COPY_N
#define TRACK_UPPER_COPY_N
#define REACTION_SIGNALING
#define RSPECIES_SIGNALING
//#define BOOST_MEM_POOL
//#define BOOL_POOL_NSIZE 65536

///Species constants
typedef unsigned short species_copy_t;
const species_copy_t max_ulim = 10000;
extern double global_time;

inline double tau() {return global_time;}

///Some constants
const double kT = 4.1; //in pN * nm
const double motorHeadGroup = 5.0; //in nm
const double linkerHeadGroup = 5.0; //in nm

///To use STL containers, libraries, etc
using namespace std;

//@{
/// Constant Reaction reactant and product numbers
/// @note - DO NOT CHANGE UNLESS YOU KNOW WHAT YOU ARE DOING!!!

/// Polymerization
#define POLYREACTANTS         2
#define POLYPRODUCTS          3

/// Depolymerization
#define DEPOLYREACTANTS       3
#define DEPOLYPRODUCTS        2

/// Linker and motor binding
#define LMBINDINGREACTANTS    3
#define LMBINDINGPRODUCTS     2

/// Linker and motor unbinding
#define LMUNBINDINGREACTANTS  2
#define LMUNBINDINGPRODUCTS   3

/// Motor walking
#define MWALKINGREACTANTS     2
#define MWALKINGPRODUCTS      2

/// Nucleation
#define NUCLEATIONREACTANTS   2
#define NUCLEATIONPRODUCTS    3

/// Destruction
#define DESTRUCTIONREACTANTS  2
#define DESTRUCTIONPRODUCTS   2

/// Aging
#define AGINGREACTANTS        1
#define AGINGPRODUCTS         1

/// Severing
#define SEVERINGREACTANTS     1
#define SEVERINGPRODUCTS      0

/// Branching
#define BRANCHINGREACTANTS    3
#define BRANCHINGPRODUCTS     2

/// Branch unbinding
#define BUNBINDINGREACTANTS   1
#define BUNBINDINGPRODUCTS    2
//@}

#endif
