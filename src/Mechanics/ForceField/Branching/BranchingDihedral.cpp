
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

#include "BranchingDihedral.h"

#include "BranchingDihedralCosine.h"

#include "BranchingPoint.h"
#include "Cylinder.h"
#include "Bead.h"

template <class BDihedralInteractionType>
void BranchingDihedral<BDihedralInteractionType>::vectorize() {
    
    beadSet = new int[n * BranchingPoint::getBranchingPoints().size()];
    kdih = new double[BranchingPoint::getBranchingPoints().size()];
    pos = new double[BranchingPoint::getBranchingPoints().size()];
    
    int i = 0;
    
    for (auto b: BranchingPoint::getBranchingPoints()) {
        
        beadSet[n * i] = b->getFirstCylinder()->getFirstBead()->_dbIndex;
        beadSet[n * i + 1] = b->getFirstCylinder()->getSecondBead()->_dbIndex;
        beadSet[n * i + 2] = b->getSecondCylinder()->getFirstBead()->_dbIndex;
        beadSet[n * i + 3] = b->getSecondCylinder()->getSecondBead()->_dbIndex;
        
        kdih[i] = b->getMBranchingPoint()->getDihedralConstant();
        pos[i] = b->getPosition();
        
        i++;
    }
}

template<class BDihedralInteractionType>
void BranchingDihedral<BDihedralInteractionType>::deallocate() {
    
    delete beadSet;
    delete kdih;
    delete pos;
}


template <class BDihedralInteractionType>
double BranchingDihedral<BDihedralInteractionType>::computeEnergy(double *coord, double *f, double d) {
    
    double U_i;
    
    if (d == 0.0)
        U_i = _FFType.energy(coord, f, beadSet, kdih, pos);
    else
        U_i = _FFType.energy(coord, f, beadSet, kdih, pos, d);
    
    return U_i;
    
}

template <class BDihedralInteractionType>
void BranchingDihedral<BDihedralInteractionType>::computeForces(double *coord, double *f) {
    
    _FFType.forces(coord, f, beadSet, kdih, pos);
}

///Template specializations
template double BranchingDihedral<BranchingDihedralCosine>::computeEnergy(double *coord, double *f, double d);
template void BranchingDihedral<BranchingDihedralCosine>::computeForces(double *coord, double *f);
template void BranchingDihedral<BranchingDihedralCosine>::vectorize();
template void BranchingDihedral<BranchingDihedralCosine>::deallocate();
