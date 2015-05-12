
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_common_h
#define M3SYM_common_h

#include <iostream>

#include "utility.h"

///Species constants
typedef unsigned short species_copy_t;
const species_copy_t max_ulim = 50000;

///Global time
extern double global_time;
inline double tau() {return global_time;}

///Some constants
const double kT = 4.1; //in pN * nm

///To use STL containers, libraries, etc
using namespace std;

/// Default empty site number in CMonomer (SpeciesBound)
/// @note - DO NOT CHANGE UNLESS YOU KNOW WHAT YOU ARE DOING!!!

#define BOUND_EMPTY       0

//@{
/// Constant Species identifiers
/// @note - DO NOT CHANGE UNLESS YOU KNOW WHAT YOU ARE DOING!!!

#define SPECIESFILAMENT   0
#define SPECIESPLUSEND    1
#define SPECIESMINUSEND   2

#define SPECIESBOUND      0
#define SPECIESLINKER     1
#define SPECIESMOTOR      2
#define SPECIESBRANCHER   3
//@}

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
