
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
            
            //do not calculate exvol for a non full length cylinder
            if(cylinder->getMCylinder()->getEqLength() !=
               SysParams::Geometry().cylinderSize) continue;
            
            for(auto &neighbor : neighborList->getNeighbors(cylinder)) {
                
                //do not calculate exvol for a non full length cylinder
                if(neighbor->getMCylinder()->getEqLength() !=
                   SysParams::Geometry().cylinderSize) continue;
                
                U_i = cylinderVolInteraction->computeEnergy(cylinder, neighbor, d);
                
                if(fabs(U_i) == numeric_limits<double>::infinity() || U_i != U_i)
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
        for(auto &cylinder : *CylinderDB::instance()) {
            
            //do not calculate exvol for a non full length cylinder
            if(cylinder->getMCylinder()->getEqLength() !=
               SysParams::Geometry().cylinderSize) continue;
            
            for(auto &neighbor : neighborList->getNeighbors(cylinder)) {
                
                //do not calculate exvol for a non full length cylinder
                if(neighbor->getMCylinder()->getEqLength() !=
                   SysParams::Geometry().cylinderSize) continue;
        
                cylinderVolInteraction->computeForces(cylinder, neighbor);
            }
        }
    }
}

void VolumeCylindricalFF::computeForcesAux() {
    
    for (auto &cylinderVolInteraction : _cylinderVolInteractionVector) {
        
        auto neighborList = cylinderVolInteraction->getNeighborList();
        for(auto &cylinder : *CylinderDB::instance()) {
            
            //do not calculate exvol for a non full length cylinder
            if(cylinder->getMCylinder()->getEqLength() !=
               SysParams::Geometry().cylinderSize) continue;
            
            for(auto &neighbor : neighborList->getNeighbors(cylinder)) {
                
                //do not calculate exvol for a non full length cylinder
                if(neighbor->getMCylinder()->getEqLength() !=
                   SysParams::Geometry().cylinderSize) continue;
                
                cylinderVolInteraction->computeForcesAux(cylinder, neighbor);
            }
        }
    }
}
