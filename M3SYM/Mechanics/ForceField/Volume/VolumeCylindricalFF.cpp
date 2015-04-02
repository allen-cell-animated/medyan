
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

#include "VolumeCylindricalFF.h"

#include "CylinderExclVolume.h"
#include "CylinderExclVolRepulsion.h"

#include "Cylinder.h"

VolumeCylindricalFF::VolumeCylindricalFF (string& type) {
    if (type == "REPULSION")
        _cylinderVolInteractionVector.emplace_back(
            new CylinderExclVolume <CylinderExclVolRepulsion>());
    else if(type == "") {}
    else {
        cout << "Volume FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
}

double VolumeCylindricalFF::computeEnergy(double d) {
    
    double U= 0;
    double U_i;
    
    for (auto &cylinderVolInteraction : _cylinderVolInteractionVector) {
        
        auto neighborList = cylinderVolInteraction->getNeighborList();
        for(auto &cylinder : *CylinderDB::instance()) {
            for(auto &cNeighbor : neighborList->getNeighbors(cylinder)) {
                
                U_i = cylinderVolInteraction->computeEnergy(cylinder, cNeighbor, d);
                
                if(U_i == numeric_limits<double>::infinity() || U != U)
                    return -1;
                else
                    U += U_i;
            }
        }
    }
    return U;
}

void VolumeCylindricalFF::computeForces() {
    
    for (auto &cylinderVolInteraction : _cylinderVolInteractionVector) {
        
        auto neighborList = cylinderVolInteraction->getNeighborList();
        for(auto &cylinder : *CylinderDB::instance())
            for(auto &cNeighbor : neighborList->getNeighbors(cylinder))
                cylinderVolInteraction->computeForces(cylinder, cNeighbor);
    }
}

void VolumeCylindricalFF::computeForcesAux() {
    
    for (auto &cylinderVolInteraction : _cylinderVolInteractionVector) {
        
        auto neighborList = cylinderVolInteraction->getNeighborList();
        for(auto &cylinder : *CylinderDB::instance())
            for(auto &cNeighbor : neighborList->getNeighbors(cylinder))
                cylinderVolInteraction->computeForcesAux(cylinder, cNeighbor);
    }
}
