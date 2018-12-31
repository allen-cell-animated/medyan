
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

#include "CylinderExclVolume.h"

#include "CylinderExclVolRepulsion.h"

#include "Cylinder.h"
#include "Bead.h"

template <class CVolumeInteractionType>
double CylinderExclVolume<CVolumeInteractionType>::computeEnergy(bool stretched) {
    
    double U = 0.0;
    double U_i=0.0;
    
    for(auto ci : Cylinder::getCylinders()) {
        
        //do not calculate exvol for a non full length cylinder
        if(!ci->isFullLength()) continue;
        
        for(auto &cn : _neighborList->getNeighbors(ci)) {
            
            //do not calculate exvol for a branching cylinder
            if(!cn->isFullLength() ||
               cn->getBranchingCylinder() == ci) continue;
            
            Bead* b1 = ci->getFirstBead();
            Bead* b2 = ci->getSecondBead();
            Bead* b3 = cn->getFirstBead();
            Bead* b4 = cn->getSecondBead();
            double kRepuls = ci->getMCylinder()->getExVolConst();
            
            U_i = _FFType.energy(b1, b2, b3, b4, kRepuls, stretched);
            
            if(fabs(U_i) == numeric_limits<double>::infinity()
               || U_i != U_i || U_i < -1.0) {
                
                //set culprits and exit
                _cylinderCulprit1 = ci;
                _cylinderCulprit2 = cn;
                
                return -1;
            }
            else
                U += U_i;
        }
    }
    
    return U;
}

template <class CVolumeInteractionType>
void CylinderExclVolume<CVolumeInteractionType>::computeForces() {

    for(auto ci : Cylinder::getCylinders()) {
        
        //do not calculate exvol for a non full length cylinder
        if(!ci->isFullLength()) continue;
        
        for(auto &cn : _neighborList->getNeighbors(ci)) {
            
            //do not calculate exvol for a branching cylinder
            if(!cn->isFullLength() ||
               cn->getBranchingCylinder() == ci) continue;
            
            Bead* b1 = ci->getFirstBead();
            Bead* b2 = ci->getSecondBead();
            Bead* b3 = cn->getFirstBead();
            Bead* b4 = cn->getSecondBead();
            double kRepuls = ci->getMCylinder()->getExVolConst();
            _FFType.forces(b1, b2, b3, b4, kRepuls);
        }
    }
}


template <class CVolumeInteractionType>
void CylinderExclVolume<CVolumeInteractionType>::computeForcesAux() {

    for(auto ci : Cylinder::getCylinders()) {
        
        //do not calculate exvol for a non full length cylinder
        if(!ci->isFullLength()) continue;
        
        for(auto &cn : _neighborList->getNeighbors(ci)) {
            //do not calculate exvol for a branching cylinder
            if(!cn->isFullLength() ||
               cn->getBranchingCylinder() == ci) continue;
            
            Bead* b1 = ci->getFirstBead();
            Bead* b2 = ci->getSecondBead();
            Bead* b3 = cn->getFirstBead();
            Bead* b4 = cn->getSecondBead();
            double kRepuls = ci->getMCylinder()->getExVolConst();
            
            _FFType.forcesAux(b1, b2, b3, b4, kRepuls);
        }
    }
}

///Template specializations
template double CylinderExclVolume<CylinderExclVolRepulsion>::computeEnergy(bool stretched);
template void CylinderExclVolume<CylinderExclVolRepulsion>::computeForces();
template void CylinderExclVolume<CylinderExclVolRepulsion>::computeForcesAux();

