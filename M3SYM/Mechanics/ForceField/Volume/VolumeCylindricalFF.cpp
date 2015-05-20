
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
    
    for (auto &interaction : _cylinderVolInteractionVector) {
        
        auto nl = interaction->getNeighborList();
        for(auto ci : Cylinder::getCylinders()) {
            
            //do not calculate exvol for a non full length cylinder
            if(ci->getMCylinder()->getEqLength() !=
               SysParams::Geometry().cylinderSize) continue;
            
            for(auto &cn : nl->getNeighbors(ci)) {
                
                //do not calculate exvol for a non full length cylinder
                if(cn->getMCylinder()->getEqLength() !=
                   SysParams::Geometry().cylinderSize) continue;
                
                U_i = interaction->computeEnergy(ci, cn, d);
                
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
    
    for (auto &interaction : _cylinderVolInteractionVector) {
        
        auto nl = interaction->getNeighborList();
        for(auto ci : Cylinder::getCylinders()) {
            
            //do not calculate exvol for a non full length cylinder
            if(!ci->isFullLength()) continue;
            
            for(auto &cn : nl->getNeighbors(ci)) {
                
                //do not calculate exvol for a non full length cylinder
                if(!cn->isFullLength()) continue;
                
                interaction->computeForces(ci, cn);
            }
        }
    }
}

void VolumeCylindricalFF::computeForcesAux() {
    
    for (auto &interaction : _cylinderVolInteractionVector) {
        
        auto nl = interaction->getNeighborList();
        for(auto ci : Cylinder::getCylinders()) {
            
            //do not calculate exvol for a non full length cylinder
            if(!ci->isFullLength()) continue;
            
            for(auto &cn : nl->getNeighbors(ci)) {
                
                //do not calculate exvol for a non full length cylinder
                if(!cn->isFullLength()) continue;
                
                interaction->computeForcesAux(ci, cn);
            }
        }
    }
}

vector<NeighborList*> VolumeCylindricalFF::getNeighborLists() {
    
    vector<NeighborList*> neighborLists;
    
    for(auto &interaction : _cylinderVolInteractionVector)
        neighborLists.push_back(interaction->getNeighborList());
    
    return neighborLists;
}

