//
//  MLinkerFF.cpp
//  Cyto
//
//  Created by Konstantin Popov on 9/5/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "MLinkerFF.h"

double LinkerFF::ComputeEnergy(double d)
{
    double U_linker = 0;
    
    for ( auto it: LinkerDB::Instance(LinkerDBKey()) )
    {
        
        for (auto linkerInteraction : _linkerInteractionVector)
        {
            U_linker += linkerInteraction->ComputeEnergy(it, d);
        }
        
    }
    
    return U_linker;
}

void LinkerFF::ComputeForces()
{
    
    for ( auto it: LinkerDB::Instance(LinkerDBKey()) )
    {
        for (auto linkerInteraction : _linkerInteractionVector)
        {
            linkerInteraction->ComputeForces(it);
        }
        
    }
    
}

void LinkerFF::ComputeForcesAux()
{
    
    for ( auto it: LinkerDB::Instance(LinkerDBKey()) )
    {
        for (auto linkerInteraction : _linkerInteractionVector)
        {
            linkerInteraction->ComputeForcesAux(it);
        }
        
    }
    
}