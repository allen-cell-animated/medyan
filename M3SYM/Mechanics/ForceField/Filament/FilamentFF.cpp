
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.0
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "FilamentFF.h"

#include "FilamentStretching.h"
#include "FilamentStretchingHarmonic.h"

#include "FilamentBending.h"
#include "FilamentBendingHarmonic.h"
#include "FilamentBendingCosine.h"

#include "Filament.h"

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

void FilamentFF::whoIsCulprit() {
    
    cout << endl;
    
    cout << "Culprit interaction = " << _culpritInteraction->getName() << endl;
    
    cout << "Printing the culprit filament..." << endl;
    _culpritInteraction->_filamentCulprit->printSelf();
    
    cout << endl;
}


double FilamentFF::computeEnergy(double d) {
    
    double U= 0;
    double U_i;
    
    for (auto &interaction : _filamentInteractionVector) {
        
        U_i = interaction->computeEnergy(d);
        
        if(U_i <= -1) {
            //set culprit and return
            _culpritInteraction = interaction.get();
            return -1;
        }
        else U += U_i;
        
    }
    return U;
}

void FilamentFF::computeForces() {
    
    for (auto &interaction : _filamentInteractionVector)
        interaction->computeForces();
}

void FilamentFF::computeForcesAux() {
    
    for (auto &interaction : _filamentInteractionVector)
        interaction->computeForcesAux();
}
