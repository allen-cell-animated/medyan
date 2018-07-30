
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

#include "BoundaryFF.h"

#include "BoundaryCylinderRepulsion.h"
#include "BoundaryCylinderRepulsionExp.h"

#include "BoundaryElement.h"
#include "Bead.h"
#include "Composite.h"

BoundaryFF::BoundaryFF (string type) {
    
    if (type == "REPULSIONEXP") {
        _boundaryInteractionVector.emplace_back(
        new BoundaryCylinderRepulsion<BoundaryCylinderRepulsionExp>());
    }
    else if(type == "") {}
    else {
        cout << "Boundary FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
}

void BoundaryFF::whoIsCulprit() {
    
    cout << endl;
    
    cout << "Culprit interaction = " << _culpritInteraction->getName() << endl;
    
    cout << "Printing the culprit boundary element..." << endl;
    _culpritInteraction->_boundaryElementCulprit->printSelf();
    
    cout << "Printing the other culprit structure..." << endl;
    _culpritInteraction->_otherCulprit->printSelf();
    
    cout << endl;
}


double BoundaryFF::computeEnergy(double d) {
    
    double U= 0;
    double U_i;
    
    for (auto &interaction : _boundaryInteractionVector) {
        
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

void BoundaryFF::computeForces() {

    for (auto &interaction : _boundaryInteractionVector)
        interaction->computeForces();
}

void BoundaryFF::computeForcesAux() {
    
    for (auto &interaction : _boundaryInteractionVector)
        interaction->computeForcesAux();
}

vector<NeighborList*> BoundaryFF::getNeighborLists() {
    
    vector<NeighborList*> neighborLists;
    
    for(auto &interaction : _boundaryInteractionVector)
        neighborLists.push_back(interaction->getNeighborList());
    
    return neighborLists;
}
