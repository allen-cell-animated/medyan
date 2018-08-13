
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "LinkerFF.h"

#include "LinkerStretching.h"
#include "LinkerStretchingHarmonic.h"

#include "Linker.h"

LinkerFF::LinkerFF (string& stretching, string& bending, string& twisting)
{
    if (stretching == "HARMONIC")
        _linkerInteractionVector.emplace_back(
        new LinkerStretching<LinkerStretchingHarmonic>());
    else if(stretching == "") {}
    else {
        cout << "Linker stretching FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
}

void LinkerFF::whoIsCulprit() {
    
    cout << endl;
    
    cout << "Culprit interaction = " << _culpritInteraction->getName() << endl;
    
    cout << "Printing the culprit linker..." << endl;
    _culpritInteraction->_linkerCulprit->printSelf();
    
    cout << endl;
}

double LinkerFF::computeEnergy(double d) {
    
    double U= 0.0;
    double U_i=0.0;
    
    for (auto &interaction : _linkerInteractionVector) {
        
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

void LinkerFF::computeForces() {
    
    for (auto &interaction : _linkerInteractionVector)
        interaction->computeForces();
}

void LinkerFF::computeForcesAux() {
    
    for (auto &interaction : _linkerInteractionVector)
        interaction->computeForcesAux();
}

