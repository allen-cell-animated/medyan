
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

std::mt19937 Rand::eng; // Default seed. Will be overwritten at initialization.
std::mt19937 Rand::engFixed(0);
std::uniform_int_distribution<int> Rand::_int_distr;
