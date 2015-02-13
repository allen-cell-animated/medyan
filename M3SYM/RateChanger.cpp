
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

#include "RateChanger.h"

float BrownianRatchet::changeRate(float bareRate, double force) {
    
    return bareRate * exp( - force * _x / kT);
}

float CatchSlipBond::changeRate(float bareRate, double force) {

    return bareRate * (_a1 * exp(- force * _x1 / kT) +
                       _a2 * exp( force * _x2 / kT));
}

float SlipBond::changeRate(float bareRate, double force) {
    
    return bareRate * exp( force * _x / kT);
}

float ExpStall::changeRate(float bareRate, double force) {
    
    return bareRate * exp( - force * _x / kT);
}

