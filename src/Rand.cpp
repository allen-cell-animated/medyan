
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
mt19937 Rand::eng(rdtsc());
#endif
uniform_int_distribution<int> Rand::_int_distr;

