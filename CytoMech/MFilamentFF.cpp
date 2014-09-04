//
//  MFilamentFF.cpp
//  Cyto
//
//  Created by Konstantin Popov on 8/19/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "MFilamentFF.h"


double FilamentFF::ComputeEnergy(double d)
{
    double U_fil = 0;
    
    for ( auto it: FilamentDB::Instance(FilamentDBKey()) )
    {
        
        for (auto filamentInteraction : _filamentInteractionVector)
        {
            U_fil += filamentInteraction->ComputeEnergy(it, d);
        }
    
    }

    return U_fil;
}

void FilamentFF::ComputeForces()
{
    
    for ( auto it: FilamentDB::Instance(FilamentDBKey()) )
    {
        
        for (auto filamentInteraction : _filamentInteractionVector)
        {
          filamentInteraction->ComputeForces(it);
        }
        
    }

}

void FilamentFF::ComputeForcesAux()
{
    
    for ( auto it: FilamentDB::Instance(FilamentDBKey()) )
    {
        
        for (auto filamentInteraction : _filamentInteractionVector)
        {
            filamentInteraction->ComputeForcesAux(it);
        }
        
    }
    
}