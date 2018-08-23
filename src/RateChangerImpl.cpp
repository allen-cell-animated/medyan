
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
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

//Qin ----------------
float BranchSlip::changeRate(float bareRate, double force) {
    
    double newRate = bareRate * exp( force * _x / kT);
    
    return newRate;
}

float MotorCatch::numBoundHeads(float onRate, float offRate,
                                double force, int numHeads) {
#ifdef PLOSFEEDBACK
    return min<double>(numHeads, numHeads * _dutyRatio + _gamma * force);
#else
    return min<double>(numHeads,numHeads * _dutyRatio + _beta * force / numHeads);
#endif
    
}

float MotorCatch::changeRate(float onRate, float offRate,
                             double numHeads, double force) {
    
    //calculate new rate
#ifdef PLOSFEEDBACK
    double k_0 = _beta * onRate /numBoundHeads(onRate, offRate, force, numHeads);

    double factor = exp(-force / (numBoundHeads(onRate, offRate, force, numHeads) * _F0));
#else
    double k_0 = onRate * (numHeads) / (exp(log((onRate + offRate) / offRate) * numHeads)
                                        - 1.0);
    
    double factor = max(0.1, exp(-force / (numBoundHeads(onRate, offRate, force, numHeads) * _F0)));
#endif
    
    double newRate = k_0 * factor;
    return newRate;
}

float MotorStall::changeRate(float onRate, float offRate,
                             double numHeads, double force) {
    
    //determine k_0
    float k_0 = ((1 - _dutyRatio) / _dutyRatio) * onRate * _stepFrac;

    //calculate new rate
#ifdef PLOSFEEDBACK
    double newRate =  max(0.0, k_0 * (_F0 - force/numHeads)
                          / (_F0 + (force / (numHeads * _alpha))));
#else
    double newRate =  max(0.0, k_0 * (_F0 - force)
                               / (_F0 + (force / (_alpha))));
#endif
    
    return newRate;
}

