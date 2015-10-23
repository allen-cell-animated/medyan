
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

#ifndef M3SYM_CGPolakRibiereMethod_h
#define M3SYM_CGPolakRibiereMethod_h

#include <cmath>
#include <numeric>
#include <algorithm>

#include "common.h"

#include "CGMethod.h"

/// The Polak-Ribiere method for conjugate gradient minimization
class PolakRibiere : public CGMethod
{
public:
    virtual void minimize(ForceFieldManager &FFM, double GRADTOL,
                          double MAXDIST, double LAMBDAMAX, bool steplimit);
};

#endif

