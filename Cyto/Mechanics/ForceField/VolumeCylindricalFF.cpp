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

VolumeCylindricalFF::VolumeCylindricalFF (std::string& Type)
{
    if (Type == "REPULSION") {_cylinderVolInteractionVector.emplace_back(new CylinderExclVolume <CylinderExclVolRepulsion>());}
}

double VolumeCylindricalFF::ComputeEnergy(double d) {
    double U_cyl= 0;
    
    for ( auto it: *CylinderDB::Instance(CylinderDBKey()) ) {
        
        for(auto &neighbor : it->getMCylinder()->getExVolNeighborsList()) {
        //neighbour list iterator to find a pair for given cylinder.
            for (auto &cylinderVolInteraction : _cylinderVolInteractionVector){
           
                U_cyl += cylinderVolInteraction->ComputeEnergy(it, neighbor->getCylinder(), d);
            }
        }
    }
    return U_cyl;
}

void VolumeCylindricalFF::ComputeForces() {
    for ( auto it: *CylinderDB::Instance(CylinderDBKey()) ) {
        
        for(auto &neighbor : it->getMCylinder()->getExVolNeighborsList()) {
            //neighbour list iterator to find a pair for given cylinder.
            for (auto &cylinderVolInteraction : _cylinderVolInteractionVector){
                
                cylinderVolInteraction->ComputeForces(it, neighbor->getCylinder());
            }
        }
    }
}

void VolumeCylindricalFF::ComputeForcesAux() {
    
    for ( auto it: *CylinderDB::Instance(CylinderDBKey()) ) {
        
        for(auto &neighbor : it->getMCylinder()->getExVolNeighborsList()) {
            //neighbour list iterator to find a pair for given cylinder.
            for (auto &cylinderVolInteraction : _cylinderVolInteractionVector){
                
                cylinderVolInteraction->ComputeForcesAux(it, neighbor->getCylinder());
            }
        }
    }
}
