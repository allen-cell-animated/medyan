
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

float LinkerCatchSlip::changeRate(float bareRate, double force) {

    return bareRate * (_a1 * exp(- force * _x1 / kT) +
                       _a2 * exp(  force * _x2 / kT));
}

float LinkerSlip::changeRate(float bareRate, double force) {
    
    return bareRate * exp( force * _x / kT);
}

float MotorPCMCatch::changeRate(float bareRate, double force) {
    
    //first determine n_b
    
    
    
    
}

float MotorHillStall::changeRate(float bareRate, double force) {
    
    return bareRate * exp( - force * _x / kT);
}

