
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
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

#include "MathFunctions.h"
using namespace mathfunc;

template <class BStretchingInteractionType>
double BranchingStretching<BStretchingInteractionType>::computeEnergy(double d) {
    
    double U = 0.0;
    double U_i=0.0;
    
    for (auto b: BranchingPoint::getBranchingPoints()) {
        
        Bead* b1 = b->getFirstCylinder()->getFirstBead();
        Bead* b2 = b->getFirstCylinder()->getSecondBead();
        Bead* b3 = b->getSecondCylinder()->getFirstBead();
        
        double kStretch = b->getMBranchingPoint()->getStretchingConstant();
        double eqLength = b->getMBranchingPoint()->getEqLength();
        double position = b->getPosition();
        
        if (d == 0.0)
            U_i = _FFType.energy(b1, b2, b3, position, kStretch, eqLength);
        else
            U_i = _FFType.energy(b1, b2, b3, position, kStretch, eqLength, d);
        
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

template <class BStretchingInteractionType>
void BranchingStretching<BStretchingInteractionType>::computeForces() {
    auto i=0;
    for (auto b: BranchingPoint::getBranchingPoints()) {
        i++;
        Bead* b1 = b->getFirstCylinder()->getFirstBead();
        Bead* b2 = b->getFirstCylinder()->getSecondBead();
        Bead* b3 = b->getSecondCylinder()->getFirstBead();

        double kStretch = b->getMBranchingPoint()->getStretchingConstant();
        double eqLength = b->getMBranchingPoint()->getEqLength();
        double position = b->getPosition();
        
        //Qin
        //_FFType.forces(b1, b2, b3, position, kStretch, eqLength);
        double f0 = _FFType.forces(b1, b2, b3, position, kStretch, eqLength);
        b->getMBranchingPoint()->stretchForce = f0;
    }
}


template <class BStretchingInteractionType>
void BranchingStretching<BStretchingInteractionType>::computeForcesAux() {
    
    for (auto b: BranchingPoint::getBranchingPoints()) {
        
        Bead* b1 = b->getFirstCylinder()->getFirstBead();
        Bead* b2 = b->getFirstCylinder()->getSecondBead();
        Bead* b3 = b->getSecondCylinder()->getFirstBead();
        
        double kStretch = b->getMBranchingPoint()->getStretchingConstant();
        double eqLength = b->getMBranchingPoint()->getEqLength();
        double position = b->getPosition();
        
        _FFType.forcesAux(b1, b2, b3, position, kStretch, eqLength);
    }
}


///Template specializations
template double
BranchingStretching<BranchingStretchingHarmonic>::computeEnergy(double d);
template void BranchingStretching<BranchingStretchingHarmonic>::computeForces();
template void BranchingStretching<BranchingStretchingHarmonic>::computeForcesAux();
