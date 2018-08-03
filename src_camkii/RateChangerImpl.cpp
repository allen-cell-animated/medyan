
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

#include <cmath>
#include <algorithm>

#include "RateChangerImpl.h"

#include "SysParams.h"

float BrownianRatchet::changeRate(float bareRate, double force) {
    
    force = min(force, 100.0);
    
    double newRate = bareRate * exp( - force * _x / kT);
    
    return newRate;
}

float CatchSlip::changeRate(float bareRate, double force) {
    
    return bareRate * (_a1 * exp(- force * _x1 / kT) +
                       _a2 * exp(  force * _x2 / kT));
}

float Slip::changeRate(float bareRate, double force) {
    
    double newRate = bareRate * exp( force * _x / kT);
    
    return newRate;
}

float LowDutyCatch::changeRate(float onRate, float offRate,
                               double numHeads, double force) {
    
    //determine N_b
    float N_b = min(double(numHeads), dutyRatio * numHeads + (force * gamma));
    
    //calculate new rate
    double newRate = beta * (offRate / N_b) * exp(-force / (N_b * _F0));
    
    return newRate;
}

float LowDutyCatchSlip::changeRate(float onRate, float offRate,
                                   double numHeads, double force) {
    
    //determine N_b
    float N_b = min(double(numHeads), dutyRatio * numHeads + (force * gamma));
    
    //calculate new rate
    double newRate = beta * (offRate / N_b) *
    (_a1 * exp(-force / (N_b * _FCatch)) +
     _a2 * exp( force / (N_b * _FSlip)));
    
    return newRate;
    
}

float LowDutyStall::changeRate(float onRate, float offRate,
                               double numHeads, double force) {
    
    //determine k_0
    float k_0 = ((1 - dutyRatio) / dutyRatio) * onRate * _stepFrac;
    
    //calculate new rate
    double newRate =  max(0.0, k_0 * (_F0 - force / numHeads)
                          / (_F0 + (force / (zeta * numHeads))));
    
    return newRate;
}