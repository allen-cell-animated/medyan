
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

#include "BubbleFF.h"

#include "BubbleCylinderRepulsion.h"
#include "BubbleCylinderRepulsionExp.h"

#include "BubbleBubbleRepulsion.h"
#include "BubbleBubbleRepulsionExp.h"

#include "Bubble.h"
#include "Bead.h"
#include "Component.h"

BubbleFF::BubbleFF (string type) {
    if (type == "REPULSIONEXP") {
        _bubbleInteractionVector.emplace_back(
        new BubbleCylinderRepulsion<BubbleCylinderRepulsionExp>());
        _bubbleInteractionVector.emplace_back(
        new BubbleBubbleRepulsion<BubbleBubbleRepulsionExp>());
    }
    else if(type == "") {}
    else {
        cout << "Bubble FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
}

void BubbleFF::whoIsCulprit() {
    
    cout << endl;
    
    cout << "Culprit interaction = " << _culpritInteraction->getName() << endl;
    
    cout << "Printing the culprit bubble..." << endl;
    _culpritInteraction->_bubbleCulprit->printSelf();
    
    cout << "Printing the other culprit structure..." << endl;
    _culpritInteraction->_otherCulprit->printSelf();
    
    cout << endl;
}


double BubbleFF::computeEnergy(double d) {
    
    double U= 0;
    double U_i;
    
    for (auto &interaction : _bubbleInteractionVector) {
        
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

void BubbleFF::computeForces() {
    
    for (auto &interaction : _bubbleInteractionVector)
        interaction->computeForces();
}

void BubbleFF::computeForcesAux() {
    
    for (auto &interaction : _bubbleInteractionVector)
        interaction->computeForcesAux();
}

vector<NeighborList*> BubbleFF::getNeighborLists() {
    
    vector<NeighborList*> neighborLists;
    
    for(auto &interaction : _bubbleInteractionVector)
        neighborLists.push_back(interaction->getNeighborList());
    
    return neighborLists;
}