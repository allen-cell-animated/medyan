
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
    
    force = min(force, 100.0);
    
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

float MotorCatch::numBoundHeads(double force, int numHeads) {
    
    return min(double(numHeads), _dutyRatio * numHeads + (force * _gamma));
    
}

float MotorCatch::changeRate(float onRate, float offRate,
                             double numBoundHeads, double force) {
    
    //calculate new rate
    double newRate = _beta * (offRate / numBoundHeads) * exp(-force / (numBoundHeads * _F0));
    
    return newRate;
}

float MotorCatchSlip::numBoundHeads(double force, int numHeads) {
    
    return min(double(numHeads), _dutyRatio * numHeads + (force * _gamma));
    
}

float MotorCatchSlip::changeRate(float onRate, float offRate,
                                 double numBoundHeads, double force) {
    
    //calculate new rate
    double newRate = _beta * (offRate / numBoundHeads) *
                     (_a1 * exp(-force / (numBoundHeads * _FCatch)) +
                      _a2 * exp( force / (numBoundHeads * _FSlip)));
    
    return newRate;
    
}

float MotorStall::changeRate(float onRate, float offRate,
                             double numHeads, double force) {
    
    //determine k_0
    float k_0 = ((1 - _dutyRatio) / _dutyRatio) * onRate * _stepFrac;
    
    //calculate new rate
    double newRate =  max(0.0, k_0 * (_F0 - force / numHeads)
                          / (_F0 + (force / (_zeta * numHeads))));
    
    return newRate;
}

