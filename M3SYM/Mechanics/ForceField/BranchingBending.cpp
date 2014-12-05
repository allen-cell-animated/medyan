//
//  BranchingBendinig.cpp
//  M3SYM
//
//  Created by Konstantin Popov on 12/3/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "BranchingBending.h"

#include "BranchingBendingCosine.h"

#include "BranchingPoint.h"
#include "Cylinder.h"
#include "Bead.h"

template <class BBendingInteractionType>
double BranchingBending<BBendingInteractionType>::computeEnergy(BranchingPoint* b, double d) {
    
    Bead* b1 = b->getFirstCylinder()->getFirstBead();
    Bead* b2 = b->getFirstCylinder()->getSecondBead();
    Bead* b3 = b->getSecondCylinder()->getFirstBead();
    Bead* b4 = b->getSecondCylinder()->getSecondBead();
    double kBend = b->getMBranchingPoint()->getBendingConstant();
    double phi0 = b->getMBranchingPoint()->getEqPhi();

    if (d == 0.0)
        return _FFType.energy(b1, b2, b3, b4, kBend, phi0);
    else
        return _FFType.energy(b1, b2, b3, b4, kBend, phi0, d);
    
}

template <class BBendingInteractionType>
void BranchingBending<BBendingInteractionType>::computeForces(BranchingPoint* b) {
   
    Bead* b1 = b->getFirstCylinder()->getFirstBead();
    Bead* b2 = b->getFirstCylinder()->getSecondBead();
    Bead* b3 = b->getSecondCylinder()->getFirstBead();
    Bead* b4 = b->getSecondCylinder()->getSecondBead();
    double kBend = b->getMBranchingPoint()->getBendingConstant();
    double phi0 = b->getMBranchingPoint()->getEqPhi();
  
    _FFType.forces(b1, b2, b3, b4, kBend, phi0);
    
}

template <class BBendingInteractionType>
void BranchingBending<BBendingInteractionType>::computeForcesAux(BranchingPoint* b) {
    
    Bead* b1 = b->getFirstCylinder()->getFirstBead();
    Bead* b2 = b->getFirstCylinder()->getSecondBead();
    Bead* b3 = b->getSecondCylinder()->getFirstBead();
    Bead* b4 = b->getSecondCylinder()->getSecondBead();
    double kBend = b->getMBranchingPoint()->getBendingConstant();
    double phi0 = b->getMBranchingPoint()->getEqPhi();
    
    _FFType.forcesAux(b1, b2, b3, b4, kBend, phi0);
}


///Template specializations
template double BranchingBending<BranchingBendingCosine>::computeEnergy(BranchingPoint* b, double d);
template void  BranchingBending<BranchingBendingCosine>::computeForces(BranchingPoint* b);
template void  BranchingBending<BranchingBendingCosine>::computeForcesAux(BranchingPoint* b);