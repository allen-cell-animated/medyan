
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

#include "CylinderVolumeFF.h"

#include "CylinderExclVolume.h"
#include "CylinderExclVolRepulsion.h"

#include "Cylinder.h"

CylinderVolumeFF::CylinderVolumeFF (string& type) {
    if (type == "REPULSION")
        _cylinderVolInteractionVector.emplace_back(
        new CylinderExclVolume <CylinderExclVolRepulsion>());
    else if(type == "") {}
    else {
        cout << "Volume FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
}

void CylinderVolumeFF::whoIsCulprit() {
    
    cout << endl;
    
    cout << "Culprit interaction = " << _culpritInteraction->getName() << endl;
    
    cout << "Printing the culprit cylinders..." << endl;
    _culpritInteraction->_cylinderCulprit1->printSelf();
    _culpritInteraction->_cylinderCulprit2->printSelf();
    
    cout << endl;
}

double CylinderVolumeFF::computeEnergy(double d) {
    
    double U= 0;
    double U_i;
    
    for (auto &interaction : _cylinderVolInteractionVector) {
        
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

void CylinderVolumeFF::computeForces() {
    
    for (auto &interaction : _cylinderVolInteractionVector)
        interaction->computeForces();
}

void CylinderVolumeFF::computeForcesAux() {
    
    for (auto &interaction : _cylinderVolInteractionVector)
        interaction->computeForcesAux();
}

vector<NeighborList*> CylinderVolumeFF::getNeighborLists() {
    
    vector<NeighborList*> neighborLists;
    
    for(auto &interaction : _cylinderVolInteractionVector)
        neighborLists.push_back(interaction->getNeighborList());
    
    return neighborLists;
}

