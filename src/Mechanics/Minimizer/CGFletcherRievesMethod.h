
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

#ifndef MEDYAN_CGFletcherRievesMethod_h
#define MEDYAN_CGFletcherRievesMethod_h

#include <cmath>
#include <numeric>
#include <algorithm>

#include "common.h"

#include "CGMethod.h"
/// The Fletcher-Rieves method for conjugate gradient minimization
    class FletcherRieves : public CGMethod {
    public:
        virtual void minimize(ForceFieldManager &FFM, floatingpoint GRADTOL,
                              floatingpoint MAXDIST, floatingpoint LAMBDAMAX, bool steplimit);
    };
#endif
