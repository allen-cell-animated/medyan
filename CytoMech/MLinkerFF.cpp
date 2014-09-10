//
//  MLinkerFF.cpp
//  Cyto
//
//  Created by Konstantin Popov on 9/5/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "MLinkerFF.h"

LinkerFF::LinkerFF (std::string Stretching, std::string Bending, std::string Twisting)
{
    if (Stretching == "HARMONIC") {_linkerInteractionVector.emplace_back(new LinkerStretching<LinkerStretchingHarmonic>());}
    //if (Bending == "HARMONIC") {_linkerInteractionVector.push_back(new LinkerBending<FilamentBendingHarmonic>());}
    //if (Twisting == "HARMONIC") {_linkerInteractionVector.push_back(new LinkerTwisting<FilamentTwistingHarmonic>());}
}

double LinkerFF::ComputeEnergy(double d)
{
    double U_linker = 0;
    
    for ( auto it: *LinkerDB::Instance(LinkerDBKey()) )
    {
        
        for (auto &linkerInteraction : _linkerInteractionVector)
        {
            U_linker += linkerInteraction.get()->ComputeEnergy(it, d);
        }
        
    }
    
    return U_linker;
}

void LinkerFF::ComputeForces()
{
    
    for ( auto it: *LinkerDB::Instance(LinkerDBKey()) )
    {
        for (auto &linkerInteraction : _linkerInteractionVector)
        {
            linkerInteraction.get()->ComputeForces(it);
        }
        
    }
    
}

void LinkerFF::ComputeForcesAux()
{
    
    for ( auto it: *LinkerDB::Instance(LinkerDBKey()) )
    {
        for (auto &linkerInteraction : _linkerInteractionVector)
        {
            linkerInteraction.get()->ComputeForcesAux(it);
        }
        
    }
    
}