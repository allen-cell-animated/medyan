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

VolumeCylindricalFF::VolumeCylindricalFF (string& type)
{
    if (type == "REPULSION") {_cylinderVolInteractionVector.emplace_back(new CylinderExclVolume <CylinderExclVolRepulsion>());}
}

double VolumeCylindricalFF::computeEnergy(double d) {
    double U_cyl= 0;
    
    for ( auto it: *CylinderDB::instance(CylinderDBKey()) ) {
        
        for(auto &neighbor : it->getMCylinder()->getExVolNeighborsList()) {
        //neighbour list iterator to find a pair for given cylinder.
            for (auto &cylinderVolInteraction : _cylinderVolInteractionVector){
                U_cyl += cylinderVolInteraction->computeEnergy(it, neighbor->getCylinder(), d);
            }
        }
    }
    return 0.5*U_cyl;
}

void VolumeCylindricalFF::computeForces() {
    
    
    for ( auto it: *CylinderDB::instance(CylinderDBKey()) ) {
        
        for(auto &neighbor : it->getMCylinder()->getExVolNeighborsList()) {
            //neighbour list iterator to find a pair for given cylinder.
            for (auto &cylinderVolInteraction : _cylinderVolInteractionVector){
                cylinderVolInteraction->computeForces(it, neighbor->getCylinder());
            }
        }
    }
}

void VolumeCylindricalFF::computeForcesAux() {
    
    for ( auto it: *CylinderDB::instance(CylinderDBKey()) ) {
        
        for(auto &neighbor : it->getMCylinder()->getExVolNeighborsList()) {
            //neighbour list iterator to find a pair for given cylinder.
            for (auto &cylinderVolInteraction : _cylinderVolInteractionVector){
                cylinderVolInteraction->computeForcesAux(it, neighbor->getCylinder());
            }
        }
    }
}
