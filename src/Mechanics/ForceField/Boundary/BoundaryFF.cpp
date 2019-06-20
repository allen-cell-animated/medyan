
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "BoundaryFF.h"

#include "BoundaryCylinderRepulsion.h"
#include "BoundaryCylinderRepulsionIn.h"
#include "BoundaryCylinderRepulsionExp.h"
#include "BoundaryCylinderRepulsionExpIn.h"

#include "BoundaryBubbleRepulsion.h"
#include "BoundaryBubbleRepulsionExp.h"

#include "BoundaryCylinderAttachment.h"
#include "BoundaryCylinderAttachmentHarmonic.h"

#include "BoundaryElement.h"
#include "Bead.h"
#include "Composite.h"

BoundaryFF::BoundaryFF (string type) {
    
    if (type == "REPULSIONEXP") {
        _boundaryInteractionVector.emplace_back(
        new BoundaryCylinderRepulsion<BoundaryCylinderRepulsionExp>());

        _boundaryInteractionVector.emplace_back(
        new BoundaryBubbleRepulsion<BoundaryBubbleRepulsionExp>());
    }
    else if(type == "REPULSIONEXPIN") {
        _boundaryInteractionVector.emplace_back(
        new BoundaryCylinderRepulsionIn<BoundaryCylinderRepulsionExpIn>());
/*        _boundaryInteractionVector.emplace_back(
        new BoundaryBubbleRepulsion<BoundaryBubbleRepulsionExp>());*/
    }
    else {
        cout << "Boundary FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
    //if pinning to boundaries
    if(SysParams::Mechanics().pinBoundaryFilaments) {
        _boundaryInteractionVector.emplace_back(
        new BoundaryCylinderAttachment<BoundaryCylinderAttachmentHarmonic>());
    }
    
    //don't change it for now
    if(SysParams::Mechanics().pinLowerBoundaryFilaments) {
        _boundaryInteractionVector.emplace_back(
        new BoundaryCylinderAttachment<BoundaryCylinderAttachmentHarmonic>());
    }
}

void BoundaryFF::vectorize() {
    
    for (auto &interaction : _boundaryInteractionVector)
        interaction->vectorize();
}

void BoundaryFF::cleanup() {
    
    for (auto &interaction : _boundaryInteractionVector)
        interaction->deallocate();
}


void BoundaryFF::whoIsCulprit() {
    
    cout << endl;
    
    cout << "Culprit interaction = " << _culpritInteraction->getName() << endl;
    
    cout << "Printing the culprit boundary element..." << endl;
    
    if(_culpritInteraction->_boundaryElementCulprit != nullptr)
        _culpritInteraction->_boundaryElementCulprit->printSelf();
    
    cout << "Printing the other culprit structure..." << endl;
    if(_culpritInteraction->_otherCulprit != nullptr)
        _culpritInteraction->_otherCulprit->printSelf();
    
    cout << endl;
}


floatingpoint BoundaryFF::computeEnergy(floatingpoint *coord, floatingpoint *f, floatingpoint d) {
    
    floatingpoint U= 0.0;
    floatingpoint U_i=0.0;
    
    for (auto &interaction : _boundaryInteractionVector) {
        
        U_i = interaction->computeEnergy(coord, f, d);
        
        if(U_i <= -1) {
            //set culprit and return
            _culpritInteraction = interaction.get();
            return -1;
        }
        else U += U_i;
        
        
    }
    
    
    return U;
}

void BoundaryFF::computeForces(floatingpoint *coord, floatingpoint *f) {

    for (auto &interaction : _boundaryInteractionVector){
        interaction->computeForces(coord, f);

    }
    
}

void BoundaryFF::computeLoadForces() {
    
    for (auto &interaction : _boundaryInteractionVector)
        interaction->computeLoadForces();
}

vector<NeighborList*> BoundaryFF::getNeighborLists() {
    
    vector<NeighborList*> neighborLists;
    
    for(auto &interaction : _boundaryInteractionVector)
        neighborLists.push_back(interaction->getNeighborList());
    
    return neighborLists;
}
