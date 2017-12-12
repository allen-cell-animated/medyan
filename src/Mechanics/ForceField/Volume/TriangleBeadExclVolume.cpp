
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2017-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "TriangleBeadExclVolume.h"

#include "TriangleBeadExclVolRepulsion.h"

#include "Triangle.h"
#include "Bead.h"

template <class TriangleBeadExclVolumeInteractionType>
double TriangleBeadExclVolume<TriangleBeadExclVolumeInteractionType>::computeEnergy(double d) {
    
    double U = 0;
    double U_i;

    // TODO: implement this
    
    for(auto ti: Triangle::getTriangles()) {
        
        for(auto &cn : _neighborList->getNeighbors(ci)) {
            
            //do not calculate exvol for a branching cylinder
            if(!cn->isFullLength() ||
               cn->getBranchingCylinder() == ci) continue;
            
            Bead* b1 = ci->getFirstBead();
            Bead* b2 = ci->getSecondBead();
            Bead* b3 = cn->getFirstBead();
            Bead* b4 = cn->getSecondBead();
            double kRepuls = ci->getMCylinder()->getExVolConst();
            
            if (d == 0.0)
                U_i = _FFType.energy(b1, b2, b3, b4, kRepuls);
            else
                U_i = _FFType.energy(b1, b2, b3, b4, kRepuls, d);
            
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

template <class TriangleBeadExclVolumeInteractionType>
void TriangleBeadExclVolume<TriangleBeadExclVolumeInteractionType>::computeForces() {

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


template <class TriangleBeadExclVolumeInteractionType>
void TriangleBeadExclVolume<TriangleBeadExclVolumeInteractionType>::computeForcesAux() {

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
template double TriangleBeadExclVolume<TriangleBeadExclVolRepulsion>::computeEnergy(double d);
template void TriangleBeadExclVolume<TriangleBeadExclVolRepulsion>::computeForces();
template void TriangleBeadExclVolume<TriangleBeadExclVolRepulsion>::computeForcesAux();

