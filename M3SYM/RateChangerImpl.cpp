
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

#include <cmath>
#include <algorithm>

#include "RateChangerImpl.h"

#include "SysParams.h"

float BrownianRatchet::changeRate(float bareRate, double force) {
    
    double newRate = bareRate * exp( - force * _x / kT);
    
    return newRate;
}

float BasicCatchSlip::changeRate(float bareRate, double force) {
    
    return bareRate * (_a1 * exp(- force * _x1 / kT) +
                       _a2 * exp(  force * _x2 / kT));
}

float BasicSlip::changeRate(float bareRate, double force) {
    
    double newRate = bareRate * exp( force * _x / kT);

    return newRate;
}

float LowDutyPCMCatch::changeRate(float onRate, float offRate,
                                  int numHeads, double force) {
    
    //determine N_b
    float N_b = min(double(numHeads), 0.1 * numHeads + (force * 0.04));

    //calculate new rate
    double rateOneSide = (offRate / N_b) * exp(-force / (N_b * _F0));
    
    double newRate = (rateOneSide * rateOneSide) / (numHeads * onRate);
    
    return newRate;

}

float LowDutyHillStall::changeRate(float onRate, float offRate,
                                   int numHeads, double force) {
    
    //determine k_0
    float k_0 = 9.0 * onRate * _stepFrac;
    
    //calculate new rate
    double newRate =  max(0.0, k_0 * (_F0 - force / numHeads)
                          / (_F0 + (force / (0.12 * numHeads))));
    
    return newRate;

    
}

