//
//  VolumeCylindricalFF.cpp
//  Cyto
//
//  Created by Konstantin Popov on 10/29/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "VolumeCylindricalFF.h"

#include "CylinderExclVolume.h"
#include "CylinderExclVolRepulsion.h"

#include "CylinderDB.h"

VolumeCylindricalFF::VolumeCylindricalFF (string& type) {
    if (type == "REPULSION") {_cylinderVolInteractionVector.emplace_back(new CylinderExclVolume <CylinderExclVolRepulsion>());}
}

double VolumeCylindricalFF::computeEnergy(double d) {
    double U_cyl= 0;
    
    for (auto &cylinderVolInteraction : _cylinderVolInteractionVector) {
        
        auto neighborList = cylinderVolInteraction->getNeighborList();
        for(auto &cylinder : *CylinderDB::instance()) {
        
            for(auto &neighbor : neighborList->getNeighbors(cylinder)) {
                Cylinder* cNeighbor = static_cast<Cylinder*>(neighbor);
                U_cyl += cylinderVolInteraction->computeEnergy(cylinder, cNeighbor, d);
            }
        }
    }
    return U_cyl;
}

void VolumeCylindricalFF::computeForces() {
    
    for (auto &cylinderVolInteraction : _cylinderVolInteractionVector) {
        
        auto neighborList = cylinderVolInteraction->getNeighborList();
        for(auto &cylinder : *CylinderDB::instance()) {
            
            for(auto &neighbor : neighborList->getNeighbors(cylinder)) {
                Cylinder* cNeighbor = static_cast<Cylinder*>(neighbor);
                cylinderVolInteraction->computeForces(cylinder, cNeighbor);
            }
        }
    }
}

void VolumeCylindricalFF::computeForcesAux() {
    
    for (auto &cylinderVolInteraction : _cylinderVolInteractionVector) {
        
        auto neighborList = cylinderVolInteraction->getNeighborList();
        for(auto &cylinder : *CylinderDB::instance()) {
            
            for(auto &neighbor : neighborList->getNeighbors(cylinder)) {
                Cylinder* cNeighbor = static_cast<Cylinder*>(neighbor);
                cylinderVolInteraction->computeForcesAux(cylinder, cNeighbor);
            }
        }
    }
}
