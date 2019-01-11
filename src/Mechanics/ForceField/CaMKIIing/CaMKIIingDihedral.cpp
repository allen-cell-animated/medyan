
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

#include "CaMKIIingDihedral.h"

#include "CaMKIIingDihedralCosine.h"

#include "CaMKIIingPoint.h"
#include "Cylinder.h"
#include "Bead.h"

template <class BDihedralInteractionType>
double CaMKIIingDihedral<BDihedralInteractionType>::computeEnergy(double d) {
    
    double U = 0;
    double U_i;
    
    for (auto b: CaMKIIingPoint::getCaMKIIingPoints()) {
        
        Bead* b1 = b->getFirstCylinder()->getFirstBead();
        Bead* b2 = b->getFirstCylinder()->getSecondBead();
        Bead* b3 = b->getSecondCylinder()->getFirstBead();
        Bead* b4 = b->getSecondCylinder()->getSecondBead();
        double kDihedr = b->getMCaMKIIingPoint()->getDihedralConstant();
        
        double position = b->getPosition();
        // TODO fix to iterate over each of the _bonds according to coordination number
#if 0
        if (d == 0.0)
            U_i = _FFType.energy(b1, b2, b3, b4, kDihedr, position);
        else
            U_i = _FFType.energy(b1, b2, b3, b4, kDihedr, position, d);
        
        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {
            
            //set culprit and return
            _camkiiingCulprit = b;
            
            return -1;
        }
        else
            U += U_i;
#endif
    }
    return U;
    
}

template <class BDihedralInteractionType>
void CaMKIIingDihedral<BDihedralInteractionType>::computeForces() {
    
    for (auto b: CaMKIIingPoint::getCaMKIIingPoints()) {
    
        Bead* b1 = b->getFirstCylinder()->getFirstBead();
        Bead* b2 = b->getFirstCylinder()->getSecondBead();
        Bead* b3 = b->getSecondCylinder()->getFirstBead();
        Bead* b4 = b->getSecondCylinder()->getSecondBead();
        double kDihedr = b->getMCaMKIIingPoint()->getDihedralConstant();
        
        double position = b->getPosition();
        // TODO fix to iterate over each of the _bonds according to coordination number
        // _FFType.forces(b1, b2, b3, b4, kDihedr, position);
    }
}

template <class BDihedralInteractionType>
void CaMKIIingDihedral<BDihedralInteractionType>::computeForcesAux() {
    
    for (auto b: CaMKIIingPoint::getCaMKIIingPoints()) {
        
        Bead* b1 = b->getFirstCylinder()->getFirstBead();
        Bead* b2 = b->getFirstCylinder()->getSecondBead();
        Bead* b3 = b->getSecondCylinder()->getFirstBead();
        Bead* b4 = b->getSecondCylinder()->getSecondBead();
        double kDihedr = b->getMCaMKIIingPoint()->getDihedralConstant();
        
        double position = b->getPosition();
        
        // TODO fix to iterate over each of the _bonds according to coordination number
        //_FFType.forcesAux(b1, b2, b3, b4, kDihedr, position);
    }
}

///Template specializations
template double CaMKIIingDihedral<CaMKIIingDihedralCosine>::computeEnergy(double d);
template void CaMKIIingDihedral<CaMKIIingDihedralCosine>::computeForces();
template void CaMKIIingDihedral<CaMKIIingDihedralCosine>::computeForcesAux();
