
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "Rand.h"
#ifdef DEBUGCONSTANTSEED
mt19937 Rand::eng(1.0);
int Rand::intcounter = 0;
int Rand::floatcounter = 0;
int Rand::chemistrycounter = 0;
#else
std::mt19937 Rand::eng; // Defined with default seed. Must be seeded at initialization.
#endif
uniform_int_distribution<int> Rand::_int_distr;
uniform_int_distribution<uint64_t> Rand::_uint64_t_distr;

