
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

#ifndef M3SYM_CGSteepestDescent_h
#define M3SYM_CGSteepestDescent_h

#include <cmath>
#include <numeric>
#include <vector>
#include <algorithm>

#include "common.h"

#include "CGMethod.h"

/// The steepest descent method for conjugate gradient minimization
class SteepestDescent : public CGMethod {
public:
    void minimize(ForceFieldManager &FFM, double GRADTOL, double MAXDIST);
};

#endif