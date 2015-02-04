
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

#include "BranchingDihedral.h"

#include "BranchingDihedralCosine.h"

#include "BranchingPoint.h"
#include "Cylinder.h"
#include "Bead.h"

template <class BDihedralInteractionType>
double BranchingDihedral<BDihedralInteractionType>::computeEnergy(
                                      BranchingPoint* b, double d) {
    
    Bead* b1 = b->getFirstCylinder()->getFirstBead();
    Bead* b2 = b->getFirstCylinder()->getSecondBead();
    Bead* b3 = b->getSecondCylinder()->getFirstBead();
    Bead* b4 = b->getSecondCylinder()->getSecondBead();
    double kDihedr = b->getMBranchingPoint()->getDihedralConstant();
    
    double position = b->getPosition();
    
    if (d == 0.0)
        return _FFType.energy(b1, b2, b3, b4, kDihedr, position);
    else
        return _FFType.energy(b1, b2, b3, b4, kDihedr, position, d);
    
}

template <class BDihedralInteractionType>
void BranchingDihedral<BDihedralInteractionType>::computeForces(BranchingPoint* b) {
    
    Bead* b1 = b->getFirstCylinder()->getFirstBead();
    Bead* b2 = b->getFirstCylinder()->getSecondBead();
    Bead* b3 = b->getSecondCylinder()->getFirstBead();
    Bead* b4 = b->getSecondCylinder()->getSecondBead();
    double kDihedr = b->getMBranchingPoint()->getDihedralConstant();
    
    double position = b->getPosition();
    
    _FFType.forces(b1, b2, b3, b4, kDihedr, position);
    
}

template <class BDihedralInteractionType>
void BranchingDihedral<BDihedralInteractionType>::computeForcesAux(BranchingPoint* b) {
    
    Bead* b1 = b->getFirstCylinder()->getFirstBead();
    Bead* b2 = b->getFirstCylinder()->getSecondBead();
    Bead* b3 = b->getSecondCylinder()->getFirstBead();
    Bead* b4 = b->getSecondCylinder()->getSecondBead();
    double kDihedr = b->getMBranchingPoint()->getDihedralConstant();
    
    double position = b->getPosition();
    
     _FFType.forcesAux(b1, b2, b3, b4, kDihedr, position);
    
}

///Template specializations
template double
BranchingDihedral<BranchingDihedralCosine>::computeEnergy(BranchingPoint* b, double d);
template void
BranchingDihedral<BranchingDihedralCosine>::computeForces(BranchingPoint* b);
template void
BranchingDihedral<BranchingDihedralCosine>::computeForcesAux(BranchingPoint* b);
