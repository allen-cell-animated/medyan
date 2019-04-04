
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

#include "CylinderVolumeFF.h"

#include "CylinderExclVolume.h"
#include "CylinderExclVolRepulsion.h"

#include "Cylinder.h"

CylinderVolumeFF::CylinderVolumeFF (string& type) {
    std::cout<<type<<endl;
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

void CylinderVolumeFF::vectorize() {
    
    for (auto &interaction : _cylinderVolInteractionVector)
        interaction->vectorize();
}

void CylinderVolumeFF::cleanup() {
    
    for (auto &interaction : _cylinderVolInteractionVector)
        interaction->deallocate();
}



floatingpoint CylinderVolumeFF::computeEnergy(floatingpoint *coord, floatingpoint *f, floatingpoint d) {

    floatingpoint U= 0.0;
    floatingpoint U_i=0.0;
    
    for (auto &interaction : _cylinderVolInteractionVector) {
        
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

#ifdef HYBRID_NLSTENCILLIST
void CylinderVolumeFF::setHNeighborLists(HybridCylinderCylinderNL* Hnl) {
    for (auto &interaction : _cylinderVolInteractionVector){
        interaction->setHNeighborList(Hnl);
    }
};
#endif

void CylinderVolumeFF::computeForces(floatingpoint *coord, floatingpoint *f) {
    
    for (auto &interaction : _cylinderVolInteractionVector)
        interaction->computeForces(coord, f);
}


vector<NeighborList*> CylinderVolumeFF::getNeighborLists() {
    
    vector<NeighborList*> neighborLists;
    
    for(auto &interaction : _cylinderVolInteractionVector)
        neighborLists.push_back(interaction->getNeighborList());
    
    return neighborLists;
}

