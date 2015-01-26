
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "BranchingPosition.h"

#include "BranchingPositionCosine.h"

#include "BranchingPoint.h"
#include "Cylinder.h"
#include "Bead.h"

template <class BPositionInteractionType>
double BranchingPosition<BPositionInteractionType>::computeEnergy(BranchingPoint* b, double d) {
    
    Bead* b1 = b->getFirstCylinder()->getFirstBead();
    Bead* b2 = b->getFirstCylinder()->getSecondBead();
    Bead* b3 = b->getSecondCylinder()->getFirstBead();
    
    double kPosition = b->getMBranchingPoint()->getPositionConstant();
    double position = b->getPosition();
    
    if (d == 0.0)
        return _FFType.energy(b1, b2, b3, position, kPosition);
    else
        return _FFType.energy(b1, b2, b3, position, kPosition, d);
    
    
}

template <class BPositionInteractionType>
void BranchingPosition<BPositionInteractionType>::computeForces(BranchingPoint* b) {
    
    Bead* b1 = b->getFirstCylinder()->getFirstBead();
    Bead* b2 = b->getFirstCylinder()->getSecondBead();
    Bead* b3 = b->getSecondCylinder()->getFirstBead();
    
    double kPosition = b->getMBranchingPoint()->getPositionConstant();
    double position = b->getPosition();
    
    _FFType.forces(b1, b2, b3, position, kPosition);
    
}


template <class BPositionInteractionType>
void BranchingPosition<BPositionInteractionType>::computeForcesAux(BranchingPoint* b) {
    
    Bead* b1 = b->getFirstCylinder()->getFirstBead();
    Bead* b2 = b->getFirstCylinder()->getSecondBead();
    Bead* b3 = b->getSecondCylinder()->getFirstBead();
    
    double kPosition = b->getMBranchingPoint()->getPositionConstant();
    double position = b->getPosition();
    
    _FFType.forcesAux(b1, b2, b3, position, kPosition);
    
}


///Template specializations
template double
BranchingPosition<BranchingPositionCosine>::computeEnergy(BranchingPoint* b, double d);
template void
BranchingPosition<BranchingPositionCosine>::computeForces(BranchingPoint* b);
template void
BranchingPosition<BranchingPositionCosine>::computeForcesAux(BranchingPoint* b);
