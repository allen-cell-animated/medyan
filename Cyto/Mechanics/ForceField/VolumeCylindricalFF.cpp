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
    //    if (Bending == "HARMONIC") {_motorGhostInteractionVector.push_back(new MotorGhostBending<MotorGhostBendingHarmonic>());}
    //    if (Twisting == "HARMONIC") {_motorGhostInteractionVector.push_back(new MotorGhostTwisting<MotorGhostTwistingHarmonic>());}
}

double VolumeCylindricalFF::ComputeEnergy(double d) {
    double U_cyl= 0;
    
    for ( auto it: *CylinderDB::Instance(CylinderDBKey()) ) {
        
        
        //neighbour list iterator to find a pair for given cylinder.
        
        for (auto &cylinderVolInteraction : _cylinderVolInteractionVector){
           
            
          //  U_cyl += cylinderVolInteraction->ComputeEnergy(it, it->neighbour, d);
        }
    }
    return U_cyl;
}

void VolumeCylindricalFF::ComputeForces() {
    for ( auto it: *CylinderDB::Instance(CylinderDBKey()) ) {
        
        
        //neighbour list iterator to find a pair for given cylinder.
        
        for (auto &cylinderVolInteraction : _cylinderVolInteractionVector){
            
            
            //  U_cyl += cylinderVolInteraction->ComputeEnergy(it, it->neighbour, d);
        }
    }
}

void VolumeCylindricalFF::ComputeForcesAux() {
    
    for ( auto it: *CylinderDB::Instance(CylinderDBKey()) ) {
        
        
        //neighbour list iterator to find a pair for given cylinder.
        
        for (auto &cylinderVolInteraction : _cylinderVolInteractionVector){
            
            
            //  U_cyl += cylinderVolInteraction->ComputeEnergy(it, it->neighbour, d);
        }
    }

}
