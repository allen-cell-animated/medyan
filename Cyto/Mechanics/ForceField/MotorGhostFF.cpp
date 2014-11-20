//
//  MotorGhostFF.cpp
//  Cyto
//
//  Created by Konstantin Popov on 9/3/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "MotorGhostFF.h"

#include "MotorGhostStretching.h"
#include "MotorGhostStretchingHarmonic.h"
#include "MotorGhostDB.h"

MotorGhostFF::MotorGhostFF (string& stretching, string& bending, string& twisting)
{
    if (stretching == "HARMONIC") {_motorGhostInteractionVector.emplace_back(new MotorGhostStretching<MotorGhostStretchingHarmonic>());}
//    if (Bending == "HARMONIC") {_motorGhostInteractionVector.push_back(new MotorGhostBending<MotorGhostBendingHarmonic>());}
//    if (Twisting == "HARMONIC") {_motorGhostInteractionVector.push_back(new MotorGhostTwisting<MotorGhostTwistingHarmonic>());}
}

double MotorGhostFF::computeEnergy(double d) {
    double U_motor = 0;
    
    for ( auto motor: *MotorGhostDB::instance())  {
        
        for (auto &motorGhostInteraction : _motorGhostInteractionVector)
            U_motor += motorGhostInteraction.get()->computeEnergy(motor, d);
    }
    return U_motor;
}

void MotorGhostFF::computeForces() {
    for (auto motor: *MotorGhostDB::instance()) {
        
        for (auto &motorGhostInteraction : _motorGhostInteractionVector)
            motorGhostInteraction.get()->computeForces(motor);
    }
}

void MotorGhostFF::computeForcesAux() {
    
    for (auto motor: *MotorGhostDB::instance()) {
        
        for (auto &motorGhostInteraction : _motorGhostInteractionVector)
            motorGhostInteraction.get()->computeForcesAux(motor);
    }
}

