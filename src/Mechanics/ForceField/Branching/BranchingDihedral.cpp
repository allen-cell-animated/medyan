
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
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
double BranchingDihedral<BDihedralInteractionType>::computeEnergy(bool stretched) {
    
    double U = 0.0;
    double U_i=0.0;
    
    for (auto b: BranchingPoint::getBranchingPoints()) {
        
        Bead* b1 = b->getFirstCylinder()->getFirstBead();
        Bead* b2 = b->getFirstCylinder()->getSecondBead();
        Bead* b3 = b->getSecondCylinder()->getFirstBead();
        Bead* b4 = b->getSecondCylinder()->getSecondBead();
        double kDihedr = b->getMBranchingPoint()->getDihedralConstant();
        
        double position = b->getPosition();
        
        U_i = _FFType.energy(b1, b2, b3, b4, kDihedr, position, stretched);
        
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

template <class BDihedralInteractionType>
void BranchingDihedral<BDihedralInteractionType>::computeForces() {
    
    for (auto b: BranchingPoint::getBranchingPoints()) {
    
        Bead* b1 = b->getFirstCylinder()->getFirstBead();
        Bead* b2 = b->getFirstCylinder()->getSecondBead();
        Bead* b3 = b->getSecondCylinder()->getFirstBead();
        Bead* b4 = b->getSecondCylinder()->getSecondBead();
        double kDihedr = b->getMBranchingPoint()->getDihedralConstant();
        
        double position = b->getPosition();
        
        //Qin
        //_FFType.forces(b1, b2, b3, b4, kDihedr, position);
        double f0 = _FFType.forces(b1, b2, b3, b4, kDihedr, position);
        b->getMBranchingPoint()->dihedralForce = f0;
        
    }
}

template <class BDihedralInteractionType>
void BranchingDihedral<BDihedralInteractionType>::computeForcesAux() {
    
    for (auto b: BranchingPoint::getBranchingPoints()) {
        
        Bead* b1 = b->getFirstCylinder()->getFirstBead();
        Bead* b2 = b->getFirstCylinder()->getSecondBead();
        Bead* b3 = b->getSecondCylinder()->getFirstBead();
        Bead* b4 = b->getSecondCylinder()->getSecondBead();
        double kDihedr = b->getMBranchingPoint()->getDihedralConstant();
        
        double position = b->getPosition();
        
        _FFType.forcesAux(b1, b2, b3, b4, kDihedr, position);
    }
}

///Template specializations
template double BranchingDihedral<BranchingDihedralCosine>::computeEnergy(bool stretched);
template void BranchingDihedral<BranchingDihedralCosine>::computeForces();
template void BranchingDihedral<BranchingDihedralCosine>::computeForcesAux();
