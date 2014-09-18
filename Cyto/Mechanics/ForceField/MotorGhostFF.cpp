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

MotorGhostFF::MotorGhostFF (std::string& Stretching, std::string& Bending, std::string& Twisting)
{
    if (Stretching == "HARMONIC") {_motorGhostInteractionVector.emplace_back(new MotorGhostStretching<MotorGhostStretchingHarmonic>());}
//    if (Bending == "HARMONIC") {_motorGhostInteractionVector.push_back(new MotorGhostBending<MotorGhostBendingHarmonic>());}
//    if (Twisting == "HARMONIC") {_motorGhostInteractionVector.push_back(new MotorGhostTwisting<MotorGhostTwistingHarmonic>());}
}

double MotorGhostFF::ComputeEnergy(double d)
{
    double U_motor = 0;
    
    for ( auto it: *MotorGhostDB::Instance(MotorGhostDBKey()) )
    {
        
        for (auto &motorGhostInteraction : _motorGhostInteractionVector)
        {
            U_motor += motorGhostInteraction.get()->ComputeEnergy(it, d);
        }
        
    }
    
    return U_motor;
}

void MotorGhostFF::ComputeForces()
{
    
    for ( auto it: *MotorGhostDB::Instance(MotorGhostDBKey()) )
    {
        
        for (auto &motorGhostInteraction : _motorGhostInteractionVector)
        {
            motorGhostInteraction.get()->ComputeForces(it);
        }
        
    }
    
}

void MotorGhostFF::ComputeForcesAux()
{
    
    for ( auto it: *MotorGhostDB::Instance(MotorGhostDBKey()) )
    {
        
        for (auto &motorGhostInteraction : _motorGhostInteractionVector)
        {
            motorGhostInteraction.get()->ComputeForcesAux(it);
        }
        
    }
    
}