
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "FilamentFF.h"

#include "FilamentDB.h"

#include "FilamentStretching.h"
#include "FilamentStretchingHarmonic.h"

#include "FilamentBending.h"
#include "FilamentBendingHarmonic.h"

FilamentFF::FilamentFF (string& stretching, string& bending, string& twisting)
{
    if (stretching == "HARMONIC")
        _filamentInteractionVector.emplace_back(new FilamentStretching<FilamentStretchingHarmonic>());
    if (bending == "HARMONIC")
        _filamentInteractionVector.emplace_back(new FilamentBending<FilamentBendingHarmonic>());
    
    //if (Twisting == "HARMONIC") {_filamentInteractionVector.push_back(new FilamentTwisting<FilamentTwistingHarmonic>());}
}


double FilamentFF::computeEnergy(double d) {
    double U_fil = 0;
    for (auto fil: *FilamentDB::instance())
        for (auto &filamentInteraction : _filamentInteractionVector)
            U_fil += filamentInteraction.get()->computeEnergy(fil, d);
    return U_fil;
}

void FilamentFF::computeForces() {
    for (auto fil: *FilamentDB::instance())
        for (auto &filamentInteraction : _filamentInteractionVector)
          filamentInteraction.get()->computeForces(fil);
}

void FilamentFF::computeForcesAux() {
    
    for (auto fil: *FilamentDB::instance())
        for (auto &filamentInteraction : _filamentInteractionVector)
            filamentInteraction.get()->computeForcesAux(fil);
}
