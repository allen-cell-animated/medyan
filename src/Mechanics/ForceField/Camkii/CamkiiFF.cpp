#ifdef CAMKII
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "CamkiiFF.h"

#include "CamkiiStretching.h"
#include "CamkiiStretchingHarmonic.h"

#include "CamkiiBending.h"
#include "CamkiiBendingCosine.h"

#include "Camkii.h"

CamkiiFF::CamkiiFF (string& stretching, string& bending) {
    
    if (stretching == "HARMONIC")
        _camkiiInteractionVector.emplace_back(
        new CamkiiStretching<CamkiiStretchingHarmonic>());
    else if(stretching == "") {}
    else {
        cout << "Camkii stretching FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
    
    if(bending == "COSINE")
        _camkiiInteractionVector.emplace_back(new CamkiiBending<CamkiiBendingCosine>());
    else if(bending == "") {}
    else {
        cout << "Camkii bending FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
}

void CamkiiFF::whoIsCulprit() {
    
    cout << endl;
    
    cout << "Culprit interaction = " << _culpritInteraction->getName() << endl;
    
    cout << "Printing the culprit camkii..." << endl;
    _culpritInteraction->_camkiiCulprit->printSelf();
    
    cout << endl;
}


double CamkiiFF::computeEnergy(double d) {
    
    double U= 0;
    double U_i;
    
    for (auto &interaction : _camkiiInteractionVector) {
        
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

void CamkiiFF::computeForces() {
    
    for (auto &interaction : _camkiiInteractionVector)
        interaction->computeForces();
}

void CamkiiFF::computeForcesAux() {
    
    for (auto &interaction : _camkiiInteractionVector)
        interaction->computeForcesAux();
}
#endif