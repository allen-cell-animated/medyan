
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.0
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include <cmath>
#include <algorithm>

#include "RateChangerImpl.h"

#include "SysParams.h"

float BrownianRatchet::changeRate(float bareRate, double force) {
    
    force = min(force, 100.0); //ceiling
    
    double newRate = bareRate * exp( - force * _x / kT);
    
    return newRate;
}

float LinkerCatchSlip::changeRate(float bareRate, double force) {
    
    return bareRate * (_a1 * exp(- force * _x1 / kT) +
                       _a2 * exp(  force * _x2 / kT));
}

float LinkerSlip::changeRate(float bareRate, double force) {
    
    double newRate = bareRate * exp( force * _x / kT);
    
    return newRate;
}

float MotorCatch::numBoundHeads(float onRate, float offRate,
                                double force, int numHeads) {
    
    return numHeads * _dutyRatio + _beta * force / numHeads;
    
}

float MotorCatch::changeRate(float onRate, float offRate,
                             double numHeads, double force) {
    
    //calculate new rate
    double k_0 = onRate * (numHeads) / (exp(log((onRate + offRate) / offRate) * numHeads) - 1);
    
    double newRate = k_0 * exp(-force / (numBoundHeads(onRate, offRate, force, numHeads) * _F0));
    return newRate;
}


float MotorStall::changeRate(float onRate, float offRate,
                             double numHeads, double force) {
    
    //determine k_0
    float k_0 = ((1 - _dutyRatio) / _dutyRatio) * onRate * _stepFrac;
    
    //calculate new rate
    double newRate =  max(0.0, k_0 * (_F0 - force)
                          / (_F0 + (force / (_alpha))));
    
    return newRate;
}
