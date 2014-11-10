//
//  FilamentFF.cpp
//  Cyto
//
//  Created by Konstantin Popov on 8/19/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "FilamentFF.h"

#include "FilamentStretching.h"
#include "FilamentStretchingHarmonic.h"
#include "FilamentBending.h"
#include "FilamentBendingHarmonic.h"
#include "FilamentDB.h"

FilamentFF::FilamentFF (string& stretching, string& bending, string& twisting)
{
    if (stretching == "HARMONIC") {_filamentInteractionVector.emplace_back(new FilamentStretching<FilamentStretchingHarmonic>());}
    if (bending == "HARMONIC") {_filamentInteractionVector.emplace_back(new FilamentBending<FilamentBendingHarmonic>());}
    //if (Twisting == "HARMONIC") {_filamentInteractionVector.push_back(new FilamentTwisting<FilamentTwistingHarmonic>());}
}


double FilamentFF::computeEnergy(double d) {
    double U_fil = 0;
    
    for ( auto it: *FilamentDB::instance(FilamentDBKey()) ) {
        
        for (auto &filamentInteraction : _filamentInteractionVector)
            U_fil += filamentInteraction.get()->computeEnergy(it, d);
    }
    return U_fil;
}

void FilamentFF::computeForces() {
    for ( auto it: *FilamentDB::instance(FilamentDBKey()) ) {
        
        for (auto &filamentInteraction : _filamentInteractionVector)
          filamentInteraction.get()->computeForces(it);
    }
}

void FilamentFF::computeForcesAux() {
    
    for ( auto it: *FilamentDB::instance(FilamentDBKey()) ) {
        
        for (auto &filamentInteraction : _filamentInteractionVector)
            filamentInteraction.get()->computeForcesAux(it);
    }
}