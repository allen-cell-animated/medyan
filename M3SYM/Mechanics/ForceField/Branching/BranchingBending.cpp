
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
    double eqTheta = b->getMBranchingPoint()->getEqTheta();

    if (d == 0.0)
        return _FFType.energy(b1, b2, b3, b4, kBend, eqTheta);
    else
        return _FFType.energy(b1, b2, b3, b4, kBend, eqTheta, d);
    
}

template <class BBendingInteractionType>
void BranchingBending<BBendingInteractionType>::computeForces(BranchingPoint* b) {
   
    Bead* b1 = b->getFirstCylinder()->getFirstBead();
    Bead* b2 = b->getFirstCylinder()->getSecondBead();
    Bead* b3 = b->getSecondCylinder()->getFirstBead();
    Bead* b4 = b->getSecondCylinder()->getSecondBead();
    double kBend = b->getMBranchingPoint()->getBendingConstant();
    double eqTheta = b->getMBranchingPoint()->getEqTheta();
  
    _FFType.forces(b1, b2, b3, b4, kBend, eqTheta);
    
}

template <class BBendingInteractionType>
void BranchingBending<BBendingInteractionType>::computeForcesAux(BranchingPoint* b) {
    
    Bead* b1 = b->getFirstCylinder()->getFirstBead();
    Bead* b2 = b->getFirstCylinder()->getSecondBead();
    Bead* b3 = b->getSecondCylinder()->getFirstBead();
    Bead* b4 = b->getSecondCylinder()->getSecondBead();
    double kBend = b->getMBranchingPoint()->getBendingConstant();
    double eqTheta = b->getMBranchingPoint()->getEqTheta();
    
    _FFType.forcesAux(b1, b2, b3, b4, kBend, eqTheta);
}


///Template specializations
template double
BranchingBending<BranchingBendingCosine>::computeEnergy(BranchingPoint* b, double d);
template void
BranchingBending<BranchingBendingCosine>::computeForces(BranchingPoint* b);
template void
BranchingBending<BranchingBendingCosine>::computeForcesAux(BranchingPoint* b);
