
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
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
mt19937 Rand::_eng(1);
#else
mt19937 Rand::_eng(rdtsc());
#endif
uniform_int_distribution<int> Rand::_int_distr;
long Rand::counter = 0;
long Rand::Dcounter = 0;
long Rand::Ncounter = 0;
