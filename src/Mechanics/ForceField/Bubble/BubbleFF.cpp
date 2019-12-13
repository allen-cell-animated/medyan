
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "BubbleFF.h"

#include "BubbleCylinderRepulsion.h"
#include "BubbleCylinderRepulsionExp.h"

#include "BubbleBubbleRepulsion.h"
#include "BubbleBubbleRepulsionExp.h"

#include "MTOCAttachment.h"
#include "MTOCAttachmentHarmonic.h"

#include "MTOCBending.h"
#include "MTOCBendingCosine.h"

#include "AFMAttachment.h"
#include "AFMAttachmentHarmonic.h"

#include "Bubble.h"
#include "Bead.h"
#include "Component.h"

BubbleFF::BubbleFF (string type, string mtoc, string afm) {
    
    //general bubble interactions
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
    
    //specifically for MTOC
    if (mtoc == "ATTACHMENTHARMONIC") {
        _bubbleInteractionVector.emplace_back(
        new MTOCAttachment<MTOCAttachmentHarmonic>());

    }
    else if(mtoc == "") {}
    else {
        cout << "MTOC FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
    
    //specifically for MTOC
    if (afm == "ATTACHMENTHARMONIC") {
        _bubbleInteractionVector.emplace_back(
        new AFMAttachment<AFMAttachmentHarmonic>());
    }
    else if(afm == "") {}
    else {
        cout << "AFM FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
}

void BubbleFF::vectorize() {

    for (auto &interaction : _bubbleInteractionVector)
    interaction->vectorize();
}

void BubbleFF::cleanup() {

    for (auto &interaction : _bubbleInteractionVector)
    interaction->deallocate();
}

void BubbleFF::whoIsCulprit() {
    
    cout << endl;
    
    cout << "Culprit interaction = " << _culpritInteraction->getName() << endl;
    
//    cout << "Printing the culprit bubble..." << endl;
//    if(_culpritInteraction->_bubbleCulprit != nullptr)
//        _culpritInteraction->_bubbleCulprit->printSelf();
    
//    cout << "Printing the other culprit structure..." << endl;
//    if(_culpritInteraction->_otherCulprit != nullptr)
//        _culpritInteraction->_otherCulprit->printSelf();
    
    cout << endl;
}


floatingpoint BubbleFF::computeEnergy(floatingpoint *coord, bool stretched) {
    
    floatingpoint U= 0.0;
    floatingpoint U_i=0.0;
    
    for (auto &interaction : _bubbleInteractionVector) {
        
        U_i = interaction->computeEnergy(coord, stretched);
        
        if(U_i <= -1) {
            //set culprit and return
            _culpritInteraction = interaction.get();
            return -1;
        }
        else U += U_i;

        #ifdef TRACKDIDNOTMINIMIZE
        if(!stretched)
            SysParams::Mininimization().tempEnergyvec.push_back(U_i);
        #endif
        
    }
    return U;
}

void BubbleFF::computeForces(floatingpoint *coord, floatingpoint *f) {
    
    for (auto &interaction : _bubbleInteractionVector)
        interaction->computeForces(coord, f);
}

//void BubbleFF::computeForcesAux() {
//
//    for (auto &interaction : _bubbleInteractionVector)
//        interaction->computeForcesAux();
//}

void BubbleFF::computeLoadForces() {
    
    for (auto &interaction : _bubbleInteractionVector)
        interaction->computeLoadForces();
}
void BubbleFF::computeLoadForce(Cylinder* c, LoadForceEnd end) const {
    for(const auto& interaction : _bubbleInteractionVector)
        interaction->computeLoadForce(c, end);
}

vector<NeighborList*> BubbleFF::getNeighborLists() {
    
    vector<NeighborList*> neighborLists;
    
    for(auto &interaction : _bubbleInteractionVector) {
        
        NeighborList* nl = interaction->getNeighborList();
        if(nl != nullptr) neighborLists.push_back(nl);
    }
    
    return neighborLists;
}

vector<string> BubbleFF::getinteractionnames(){
    vector<string> temp;
    for (auto &interaction : _bubbleInteractionVector) {
        temp.push_back(interaction->getName());
    }
    return temp;
}

