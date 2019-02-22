
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_common_h
#define MEDYAN_common_h

#include <iostream>
#include <vector>
#include <boost/signals2/shared_connection_block.hpp>

#include "utility.h"

#ifdef CUDAACCL
#include <cuda.h>
#include <cuda_runtime.h>
#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#else
#define CUDA_HOSTDEV
#endif

#endif

///Species constants
typedef unsigned int species_copy_t;
const species_copy_t max_ulim = 1000000;

///Global time
extern double global_time;

inline double tau() {return global_time;}
inline void resetglobaltime() {global_time=0.0;}
///Some constants
const double kT = 4.1; //in pN * nm

const int cylinder_cache = 500;
const int bead_cache = 1000;//total number of beads that can be appended before
// revectorization

///To use STL containers, libraries, etc
using namespace std;

///Boost typedef
typedef boost::signals2::shared_connection_block ConnectionBlock;

///Num filament types maximum
#define MAX_FILAMENT_TYPES 10

//@{
/// Constant Species index identifiers
/// @note - DO NOT CHANGE!!!

#define SPECIESFILAMENT       0
#define SPECIESPLUSEND        1
#define SPECIESMINUSEND       2

#define SPECIESBOUND          0
#define SPECIESLINKER         1
#define SPECIESMOTOR          2
#define SPECIESBRANCHER       3
//@}

//@{
/// Constant Reaction reactant and product numbers
/// @note - DO NOT CHANGE!!!

/// Polymerization
#define POLYREACTANTS         2
#define POLYPRODUCTS          2

/// Depolymerization
#define DEPOLYREACTANTS       2
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
