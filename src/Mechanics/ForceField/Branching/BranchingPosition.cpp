
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "BranchingPosition.h"

#include "BranchingPositionCosine.h"

#include "BranchingPoint.h"
#include "Cylinder.h"
#include "Bead.h"

template <class BPositionInteractionType>
double BranchingPosition<BPositionInteractionType>::computeEnergy(double d) {
    
    double U = 0.0;
    double U_i=0.0;
    
    for (auto b: BranchingPoint::getBranchingPoints()) {
        
        Bead* b1 = b->getFirstCylinder()->getFirstBead();
        Bead* b2 = b->getFirstCylinder()->getSecondBead();
        Bead* b3 = b->getSecondCylinder()->getFirstBead();
        
        double kPosition = b->getMBranchingPoint()->getPositionConstant();
        double position = b->getPosition();
        
        if (d == 0.0)
            U_i = _FFType.energy(b1, b2, b3, kPosition, position);
        else
            U_i = _FFType.energy(b1, b2, b3, kPosition, position, d);
        
        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {
            
            //set culprit and return
            _branchingCulprit = b;
            
            return -1;
        }
        else
            U += U_i;
    }
    
    return U;
}

template <class BPositionInteractionType>
void BranchingPosition<BPositionInteractionType>::computeForces() {
    
    for (auto b: BranchingPoint::getBranchingPoints()) {
    
        Bead* b1 = b->getFirstCylinder()->getFirstBead();
        Bead* b2 = b->getFirstCylinder()->getSecondBead();
        Bead* b3 = b->getSecondCylinder()->getFirstBead();
        
        double kPosition = b->getMBranchingPoint()->getPositionConstant();
        double position = b->getPosition();
        
        _FFType.forces(b1, b2, b3, kPosition, position);
    }
    
}


template <class BPositionInteractionType>
void BranchingPosition<BPositionInteractionType>::computeForcesAux() {
    
    for (auto b: BranchingPoint::getBranchingPoints()) {
        
        Bead* b1 = b->getFirstCylinder()->getFirstBead();
        Bead* b2 = b->getFirstCylinder()->getSecondBead();
        Bead* b3 = b->getSecondCylinder()->getFirstBead();
        
        double kPosition = b->getMBranchingPoint()->getPositionConstant();
        double position = b->getPosition();
        
        _FFType.forcesAux(b1, b2, b3, kPosition, position);
    }
}


///Template specializations
template double BranchingPosition<BranchingPositionCosine>::computeEnergy(double d);
template void BranchingPosition<BranchingPositionCosine>::computeForces();
template void BranchingPosition<BranchingPositionCosine>::computeForcesAux();
