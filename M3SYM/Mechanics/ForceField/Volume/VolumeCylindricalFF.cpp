
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

void VolumeCylindricalFF::whoIsCulprit() {
    
    cout << endl;
    
    cout << "Printing the culprit cylinders..." << endl;
    
    _cylinderCulprit1->printInfo();
    _cylinderCulprit2->printInfo();
    
    cout << endl;
}

double VolumeCylindricalFF::computeEnergy(double d) {
    
    double U= 0;
    double U_i;
    
    for (auto &interaction : _cylinderVolInteractionVector) {
        
        auto nl = interaction->getNeighborList();
        for(auto ci : Cylinder::getCylinders()) {
            
            //do not calculate exvol for a non full length cylinder
            if(!ci->isFullLength()) continue;
            
            for(auto &cn : nl->getNeighbors(ci)) {
                
                //do not calculate exvol for a branching cylinder
                if(cn->getBranchingCylinder() == ci) continue;
                
                U_i = interaction->computeEnergy(ci, cn, d);
                
                if(fabs(U_i) == numeric_limits<double>::infinity() || U_i != U_i) {
                    
                    //set culprits and exit
                    _cylinderCulprit1 = ci;
                    _cylinderCulprit2 = cn;
                    
                    return -1;
                }
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
                
                //do not calculate exvol for a branching cylinder
                if(cn->getBranchingCylinder() == ci) continue;
                
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
                //and for a branching cylinder
                if(!cn->isFullLength() ||
                   ci->getBranchingCylinder() == cn) continue;
                
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

