
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "FilamentFF.h"

#include "FilamentDB.h"

#include "FilamentStretching.h"
#include "FilamentStretchingHarmonic.h"

#include "FilamentBending.h"
#include "FilamentBendingHarmonic.h"
#include "FilamentBendingCosine.h"

FilamentFF::FilamentFF (string& stretching, string& bending, string& twisting) {
    
    if (stretching == "HARMONIC")
        _filamentInteractionVector.emplace_back(
            new FilamentStretching<FilamentStretchingHarmonic>());
    else if(stretching == "") {}
    else {
        cout << "Filament stretching FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
    
    if (bending == "HARMONIC")
        _filamentInteractionVector.emplace_back(
            new FilamentBending<FilamentBendingHarmonic>());
    else if(bending == "COSINE")
        _filamentInteractionVector.emplace_back(
            new FilamentBending<FilamentBendingCosine>());
    else if(bending == "") {}
    else {
        cout << "Filament bending FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
}


double FilamentFF::computeEnergy(double d) {
    
    double U = 0;
    double U_i;
    
    for (auto fil: *FilamentDB::instance()) {
        for (auto &filamentInteraction : _filamentInteractionVector) {
    
            U_i = filamentInteraction.get()->computeEnergy(fil, d);
            
            if(fabs(U_i) == numeric_limits<double>::infinity() || U_i != U_i)
                return -1;
            else
                U += U_i;
        }
    }
    return U;
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
