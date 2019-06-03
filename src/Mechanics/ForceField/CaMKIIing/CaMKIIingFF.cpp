
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

#include "CaMKIIingFF.h"

#include "CaMKIIingStretching.h"
#include "CaMKIIingStretchingHarmonic.h"

#include "CaMKIIingBending.h"
#include "CaMKIIingBendingCosine.h"

#include "CaMKIIingDihedral.h"
#include "CaMKIIingDihedralCosine.h"

#include "CaMKIIingPosition.h"
#include "CaMKIIingPositionCosine.h"

#include "CaMKIIingPoint.h"

CaMKIIingFF::CaMKIIingFF(string& stretching, string& bending,
                         string& dihedral, string& position)
{
    if(stretching == "HARMONIC")
        _camkiiingInteractionVector.emplace_back(
        new CaMKIIingStretching<CaMKIIingStretchingHarmonic>());
    else if(stretching == "") {}
    else {
        cout << "CaMKIIing stretching FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
    
    if(bending == "COSINE")
        _camkiiingInteractionVector.emplace_back(
        new CaMKIIingBending<CaMKIIingBendingCosine>());
    else if(bending == "") {}
    else {
        cout << "CaMKIIing bending FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
    
    if(dihedral == "COSINE")
      _camkiiingInteractionVector.emplace_back(
      new CaMKIIingDihedral<CaMKIIingDihedralCosine>());
    else if(dihedral == "") {}
    else {
        cout << "CaMKIIing dihedral FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
    
    if(position == "COSINE")
      _camkiiingInteractionVector.emplace_back(new
          CaMKIIingPosition<CaMKIIingPositionCosine>());
    else if(position == "") {}
    else {
        cout << "CaMKIIing position FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }

}

void CaMKIIingFF::whoIsCulprit() {
    
    cout << endl;
    
    cout << "Culprit interaction = " << _culpritInteraction->getName() << endl;
    
    cout << "Printing the culprit camkiiing point..." << endl;
    _culpritInteraction->_camkiiingCulprit->printSelf();
    
    cout << endl;
}


double CaMKIIingFF::computeEnergy(double d) {
    
    double U= 0;
    double U_i;
    
    for (auto &interaction : _camkiiingInteractionVector) {
        
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

void CaMKIIingFF::computeForces() {
    
    for (auto &interaction : _camkiiingInteractionVector) {
        interaction->computeForces();
    }
}

void CaMKIIingFF::computeForcesAux() {
    
    for (auto &interaction : _camkiiingInteractionVector)
        interaction->computeForcesAux();
}
