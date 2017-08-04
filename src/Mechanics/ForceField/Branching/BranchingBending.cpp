
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

#include "BranchingBending.h"

#include "BranchingBendingCosine.h"

#include "BranchingPoint.h"
#include "Cylinder.h"
#include "Bead.h"

template <class BBendingInteractionType>
void BranchingBending<BBendingInteractionType>::vectorize() {
    
    beadSet = new int[n * BranchingPoint::getBranchingPoints().size()];
    kbend = new double[BranchingPoint::getBranchingPoints().size()];
    eqt = new double[BranchingPoint::getBranchingPoints().size()];
    
    int i = 0;
    
    for (auto b: BranchingPoint::getBranchingPoints()) {
        
        beadSet[n * i] = b->getFirstCylinder()->getFirstBead()->_dbIndex;
        beadSet[n * i + 1] = b->getFirstCylinder()->getSecondBead()->_dbIndex;
        beadSet[n * i + 2] = b->getSecondCylinder()->getFirstBead()->_dbIndex;
        
        kbend[i] = b->getMBranchingPoint()->getStretchingConstant();
        eqt[i] = b->getMBranchingPoint()->getEqTheta();
        
        i++;
    }
}

template<class BBendingInteractionType>
void BranchingBending<BBendingInteractionType>::deallocate() {
    
    delete beadSet;
    delete kbend;
    delete eqt;
}



template <class BBendingInteractionType>
double BranchingBending<BBendingInteractionType>::computeEnergy(double *coord, double *f, double d) {
    
    double U_i;
    
    if (d == 0.0)
        U_i = _FFType.energy(coord, f, beadSet, kbend, eqt);
    else
        U_i = _FFType.energy(coord, f, beadSet, kbend, eqt, d);
    
    return U_i;
}

template <class BBendingInteractionType>
void BranchingBending<BBendingInteractionType>::computeForces(double *coord, double *f) {
   
    _FFType.forces(coord, f, beadSet, kbend, eqt);
}

///Template specializations
template double BranchingBending<BranchingBendingCosine>::computeEnergy(double *coord, double *f, double d);
template void BranchingBending<BranchingBendingCosine>::computeForces(double *coord, double *f);
