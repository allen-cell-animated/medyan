
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

#include "ForceFieldManager.h"
#include <algorithm>

#include "SubSystem.h"

void ForceFieldManager::updateGeometries(bool calcDerivative, double d) {
    for(auto g : _subSystem->getGeometrics()) g->updateGeometry(calcDerivative, d);
}

void ForceFieldManager::computeForces() {
    
    //reset
    resetForces();
    
    //recompute
    for(auto &f : _forceFields) f->computeForces();
    
    //copy to auxs
    for(auto b: Bead::getBeads())
        b->forceAux = b->forceAuxP = b->force;
}

void ForceFieldManager::computeForcesAux() {
    
    //reset just aux
    resetForcesAux();
    
    //recompute
    for(auto &f : _forceFields)
        f->computeForcesAux();
}

void ForceFieldManager::computeForcesAuxP() {
    
    //copy to auxp
    for(auto b: Bead::getBeads())
        b->forceAuxP = b->forceAux;
}

void ForceFieldManager::computeLoadForces() {
    
    for(auto &f : _forceFields)
        f->computeLoadForces();
    
    //reset lfip and lfim as well
    for(auto b: Bead::getBeads()) {
        b->lfip = 0;
        b->lfim = 0;
    }
}


void ForceFieldManager::resetForces() {
    
    for(auto b: Bead::getBeads()) {
        b->force.assign (3, 0); //Set force to zero;
        std::fill(b->loadForcesP.begin(), b->loadForcesP.end(), 0.0); //Set load force to zero
        std::fill(b->loadForcesM.begin(), b->loadForcesM.end(), 0.0); //Set load force to zero
//        std::memset((void*)(&b->loadForcesP[0]), 0, sizeof(double)*(b->loadForcesP).size());  //Set load force to zero;
//        std::memset((void*)(&b->loadForcesM[0]), 0, sizeof(double)*(b->loadForcesM).size());  //Set load force to zero;
    }
}

void ForceFieldManager::resetForcesAux() {
    
    for(auto b: Bead::getBeads())
        b->forceAux.assign (3, 0); //Set forceAux to zero;
}
