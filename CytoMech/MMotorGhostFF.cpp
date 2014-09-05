//
//  MMotorGhostFF.cpp
//  Cyto
//
//  Created by Konstantin Popov on 9/3/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "MMotorGhostFF.h"


double MotorGhostFF::ComputeEnergy(double d)
{
    double U_motor = 0;
    
    for ( auto it: MotorGhostDB::Instance(MotorghostDBKey()) )
    {
        
        for (auto motorGhostInteraction : _motorGhostInteractionVector)
        {
            U_motor += motorGhostInteraction->ComputeEnergy(it, d);
        }
        
    }
    
    return U_motor;
}

void MotorGhostFF::ComputeForces()
{
    
    for ( auto it: MotorGhostDB::Instance(MotorghostDBKey()) )
    {
        
        for (auto motorGhostInteraction : _motorGhostInteractionVector)
        {
            motorGhostInteraction->ComputeForces(it);
        }
        
    }
    
}

void MotorGhostFF::ComputeForcesAux()
{
    
    for ( auto it: MotorGhostDB::Instance(MotorghostDBKey()) )
    {
        
        for (auto motorGhostInteraction : _motorGhostInteractionVector)
        {
            motorGhostInteraction->ComputeForcesAux(it);
        }
        
    }
    
}