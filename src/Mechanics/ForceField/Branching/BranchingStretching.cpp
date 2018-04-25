
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

#include "MathFunctions.h"
using namespace mathfunc;

template <class BStretchingInteractionType>
double BranchingStretching<BStretchingInteractionType>::computeEnergy(double d) {
    
    double U = 0;
    double U_i;
    
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
        Bead* b4 = b->getSecondCylinder()->getSecondBead();
//        std::cout<<i<<" "<<b->getFirstCylinder()->getID()<<" "<<twoPointDistance(b1->coordinate, b2->coordinate)<<" "<<b->getSecondCylinder()->getID()<<" "<<twoPointDistance(b3->coordinate, b4->coordinate)<<endl;

        double kStretch = b->getMBranchingPoint()->getStretchingConstant();
        double eqLength = b->getMBranchingPoint()->getEqLength();
        double position = b->getPosition();
        
        _FFType.forces(b1, b2, b3, position, kStretch, eqLength);
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
