
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

#ifndef M3SYM_CGFletcherRievesMethod_h
#define M3SYM_CGFletcherRievesMethod_h

#include <cmath>
#include <numeric>
#include <vector>
#include <algorithm>

#include "common.h"

#include "CGMethod.h"

/// The Fletcher-Rieves method for conjugate gradient minimization
class FletcherRieves : public CGMethod {
public:
   void minimize(ForceFieldManager &FFM, double GRADTOL);
};

#endif
