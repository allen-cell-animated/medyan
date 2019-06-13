
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

#include "Rand.h"


#ifdef RAND_STATIC_SEED
default_random_engine Rand::_eng(RAND_STATIC_SEED);
uniform_int_distribution<int> Rand::_int_distr;
#else
mt19937 Rand::_eng(rdtsc());
uniform_int_distribution<int> Rand::_int_distr;
#endif
