
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
double BranchingBending<BBendingInteractionType>::computeEnergy(double d) {
    
    double U = 0;
    double U_i;
    
    for (auto b: BranchingPoint::getBranchingPoints()) {
        
        Bead* b1 = b->getFirstCylinder()->getFirstBead();
        Bead* b2 = b->getFirstCylinder()->getSecondBead();
        Bead* b3 = b->getSecondCylinder()->getFirstBead();
        Bead* b4 = b->getSecondCylinder()->getSecondBead();
        double kBend = b->getMBranchingPoint()->getBendingConstant();
        double eqTheta = b->getMBranchingPoint()->getEqTheta();
        
        if (d == 0.0)
            U_i = _FFType.energy(b1, b2, b3, b4, kBend, eqTheta);
        else
            U_i = _FFType.energy(b1, b2, b3, b4, kBend, eqTheta, d);
        
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

template <class BBendingInteractionType>
void BranchingBending<BBendingInteractionType>::computeForces() {
   
    for (auto b: BranchingPoint::getBranchingPoints()) {
    
        Bead* b1 = b->getFirstCylinder()->getFirstBead();
        Bead* b2 = b->getFirstCylinder()->getSecondBead();
        Bead* b3 = b->getSecondCylinder()->getFirstBead();
        Bead* b4 = b->getSecondCylinder()->getSecondBead();
        double kBend = b->getMBranchingPoint()->getBendingConstant();
        double eqTheta = b->getMBranchingPoint()->getEqTheta();
      
        _FFType.forces(b1, b2, b3, b4, kBend, eqTheta);
    }
}

template <class BBendingInteractionType>
void BranchingBending<BBendingInteractionType>::computeForcesAux() {
    
    for (auto b: BranchingPoint::getBranchingPoints()) {
        
        Bead* b1 = b->getFirstCylinder()->getFirstBead();
        Bead* b2 = b->getFirstCylinder()->getSecondBead();
        Bead* b3 = b->getSecondCylinder()->getFirstBead();
        Bead* b4 = b->getSecondCylinder()->getSecondBead();
        double kBend = b->getMBranchingPoint()->getBendingConstant();
        double eqTheta = b->getMBranchingPoint()->getEqTheta();
        
        _FFType.forcesAux(b1, b2, b3, b4, kBend, eqTheta);
    }
}


///Template specializations
template double BranchingBending<BranchingBendingCosine>::computeEnergy(double d);
template void BranchingBending<BranchingBendingCosine>::computeForces();
template void BranchingBending<BranchingBendingCosine>::computeForcesAux();
