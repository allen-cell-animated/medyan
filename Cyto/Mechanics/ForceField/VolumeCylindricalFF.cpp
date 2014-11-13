//
//  VolumeCylindricalFF.cpp
//  Cyto
//
//  Created by Konstantin Popov on 10/29/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "VolumeCylindricalFF.h"

#include "CylinderDB.h"

#include "CylinderExclVolume.h"
#include "CylinderExclVolRepulsion.h"


VolumeCylindricalFF::VolumeCylindricalFF (string& type) {
    if (type == "REPULSION") {_cylinderVolInteractionVector.emplace_back(new CylinderExclVolume <CylinderExclVolRepulsion>());}
}

double VolumeCylindricalFF::computeEnergy(double d) {
    double U_cyl= 0;
    
    for (auto &cylinderVolInteraction : _cylinderVolInteractionVector) {
        
        auto neighborList = cylinderVolInteraction->getNeighborList();
        for(auto &cylinder : *CylinderDB::instance(CylinderDBKey())) {
        
            for(auto &neighbor : neighborList->getNeighbors(cylinder)) {
                Cylinder* c1 = static_cast<Cylinder*>(cylinder);
                Cylinder* c2 = static_cast<Cylinder*>(neighbor);
                U_cyl += cylinderVolInteraction->computeEnergy(c1, c2, d);
            }
        }
    }
    return U_cyl;
}

void VolumeCylindricalFF::computeForces() {
    
    for (auto &cylinderVolInteraction : _cylinderVolInteractionVector) {
        
        auto neighborList = cylinderVolInteraction->getNeighborList();
        for(auto &cylinder : *CylinderDB::instance(CylinderDBKey())) {
            
            for(auto &neighbor : neighborList->getNeighbors(cylinder)) {
                Cylinder* c1 = static_cast<Cylinder*>(cylinder);
                Cylinder* c2 = static_cast<Cylinder*>(neighbor);
                cylinderVolInteraction->computeForces(c1, c2);
            }
        }
    }
}

void VolumeCylindricalFF::computeForcesAux() {
    
    for (auto &cylinderVolInteraction : _cylinderVolInteractionVector) {
        
        auto neighborList = cylinderVolInteraction->getNeighborList();
        for(auto &cylinder : *CylinderDB::instance(CylinderDBKey())) {
            
            for(auto &neighbor : neighborList->getNeighbors(cylinder)) {
                Cylinder* c1 = static_cast<Cylinder*>(cylinder);
                Cylinder* c2 = static_cast<Cylinder*>(neighbor);
                cylinderVolInteraction->computeForcesAux(c1, c2);
            }
        }
    }
}
