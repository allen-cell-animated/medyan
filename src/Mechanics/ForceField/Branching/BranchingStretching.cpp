
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

#include "BranchingStretching.h"

#include "BranchingStretchingHarmonic.h"

#include "BranchingPoint.h"
#include "Cylinder.h"
#include "Bead.h"


template <class BStretchingInteractionType>
void BranchingStretching<BStretchingInteractionType>::vectorize() {
    
    beadSet = new int[n * BranchingPoint::getBranchingPoints().size()];
    kstr = new double[BranchingPoint::getBranchingPoints().size()];
    eql = new double[BranchingPoint::getBranchingPoints().size()];
    pos = new double[BranchingPoint::getBranchingPoints().size()];
    
    
    int i = 0;
    
    for (auto b: BranchingPoint::getBranchingPoints()) {
        
        beadSet[n * i] = b->getFirstCylinder()->getFirstBead()->_dbIndex;
        beadSet[n * i + 1] = b->getFirstCylinder()->getSecondBead()->_dbIndex;
        beadSet[n * i + 2] = b->getSecondCylinder()->getFirstBead()->_dbIndex;
        
        kstr[i] = b->getMBranchingPoint()->getStretchingConstant();
        eql[i] = b->getMBranchingPoint()->getEqLength();
        pos[i] = b->getPosition();
        
        i++;
    }
}

template<class BStretchingInteractionType>
void BranchingStretching<BStretchingInteractionType>::deallocate() {
    
    delete beadSet;
    delete kstr;
    delete eql;
    delete pos;
}


template <class BStretchingInteractionType>
double BranchingStretching<BStretchingInteractionType>::computeEnergy(double *coord, double *f, double d) {
    
    double U_i;
    
    if (d == 0.0)
        U_i = _FFType.energy(coord, f, beadSet, kstr, eql, pos);
    else
        U_i = _FFType.energy(coord, f, beadSet, kstr, eql, pos, d);
    
    return U_i;
}

template <class BStretchingInteractionType>
void BranchingStretching<BStretchingInteractionType>::computeForces(double *coord, double *f) {
    
    _FFType.forces(coord, f, beadSet, kstr, eql, pos);
}



///Template specializations
template double
BranchingStretching<BranchingStretchingHarmonic>::computeEnergy(double *coord, double *f, double d);
template void BranchingStretching<BranchingStretchingHarmonic>::computeForces(double *coord, double *f);
