
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
void BranchingPosition<BPositionInteractionType>::vectorize() {
    
    beadSet = new int[n * BranchingPoint::getBranchingPoints().size()];
    kpos = new double[BranchingPoint::getBranchingPoints().size()];
    pos = new double[BranchingPoint::getBranchingPoints().size()];

    int i = 0;
    
    for (auto b: BranchingPoint::getBranchingPoints()) {
        
        beadSet[n * i] = b->getFirstCylinder()->getFirstBead()->_dbIndex;
        beadSet[n * i + 1] = b->getFirstCylinder()->getSecondBead()->_dbIndex;
        beadSet[n * i + 2] = b->getSecondCylinder()->getFirstBead()->_dbIndex;
        
        kpos[i] = b->getMBranchingPoint()->getPositionConstant();
        pos[i] = b->getPosition();
        
        i++;
    }
}

template<class BPositionInteractionType>
void BranchingPosition<BPositionInteractionType>::deallocate() {
    
    delete beadSet;
    delete kpos;
    delete pos;
}

template <class BPositionInteractionType>
double BranchingPosition<BPositionInteractionType>::computeEnergy(double *coord, double *f, double d) {
    
    double U_i;
    
    if (d == 0.0)
        U_i = _FFType.energy(coord, f, beadSet, kpos, pos);
    else
        U_i = _FFType.energy(coord, f, beadSet, kpos, pos, d);
    
    return U_i;
}

template <class BPositionInteractionType>
void BranchingPosition<BPositionInteractionType>::computeForces(double *coord, double *f) {
    
    _FFType.energy(coord, f, beadSet, kpos, pos);
    
}


///Template specializations
template double BranchingPosition<BranchingPositionCosine>::computeEnergy(double *coord, double *f, double d);
template void BranchingPosition<BranchingPositionCosine>::computeForces(double *coord, double *f);
template void BranchingPosition<BranchingPositionCosine>::vectorize();
template void BranchingPosition<BranchingPositionCosine>::deallocate();
