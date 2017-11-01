
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

float MotorCatch::changeRate(float onRate, float offRate, double force) {
    
    //calculate new rate
    double factor = 0.9 * exp(-force / _F0_catch) + 0.1 * exp(force/ _F0_slip);
    
    double newRate = _k_0 * factor;
    return newRate;
}
float MotorStall::changeRate(float onRate, float offRate, double force) {

    
    //determine k_0 (FOR MYOSIN-ISOFORMS)
    float k_0 = v_0 / _stepFrac;
    
    
    //calculate new rate
    double newRate =  max(0.0, k_0 * (1 - force / _F0));
    
    return newRate;
}
