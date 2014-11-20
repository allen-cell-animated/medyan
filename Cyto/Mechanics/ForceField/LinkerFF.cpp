//
//  LinkerFF.cpp
//  Cyto
//
//  Created by Konstantin Popov on 9/5/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "LinkerFF.h"

#include "LinkerStretching.h"
#include "LinkerStretchingHarmonic.h"
#include "LinkerDB.h"

LinkerFF::LinkerFF (string& stretching, string& bending, string& twisting)
{
    if (stretching == "HARMONIC") {_linkerInteractionVector.emplace_back(new LinkerStretching<LinkerStretchingHarmonic>());}
    //if (Bending == "HARMONIC") {_linkerInteractionVector.push_back(new LinkerBending<FilamentBendingHarmonic>());}
    //if (Twisting == "HARMONIC") {_linkerInteractionVector.push_back(new LinkerTwisting<FilamentTwistingHarmonic>());}
}

double LinkerFF::computeEnergy(double d) {
    double U_linker = 0;
    
    for (auto linker: *LinkerDB::instance()) {
        
        for (auto &linkerInteraction : _linkerInteractionVector)
            U_linker += linkerInteraction.get()->computeEnergy(linker, d);
    }
    return U_linker;
}

void LinkerFF::computeForces() {
    
    for (auto linker: *LinkerDB::instance()) {
        
        for (auto &linkerInteraction : _linkerInteractionVector)
            linkerInteraction.get()->computeForces(linker);
    }
}

void LinkerFF::computeForcesAux() {
    
    for (auto linker: *LinkerDB::instance()) {
    
        for (auto &linkerInteraction : _linkerInteractionVector)
            linkerInteraction.get()->computeForcesAux(linker);
    }
}

