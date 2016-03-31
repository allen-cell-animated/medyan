
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef M3SYM_CGSteepestDescent_h
#define M3SYM_CGSteepestDescent_h

#include <cmath>
#include <numeric>
#include <algorithm>

#include "common.h"

#include "CGMethod.h"

/// The steepest descent method for conjugate gradient minimization
class SteepestDescent : public CGMethod {
public:
    virtual void minimize(ForceFieldManager &FFM, double GRADTOL,
                          double MAXDIST, double LAMBDAMAX, bool steplimit);
};

#endif