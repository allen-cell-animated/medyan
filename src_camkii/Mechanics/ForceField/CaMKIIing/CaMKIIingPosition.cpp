
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

#include "CaMKIIingPosition.h"

#include "CaMKIIingPositionCosine.h"

#include "CaMKIIingPoint.h"
#include "Cylinder.h"
#include "Bead.h"

template <class BPositionInteractionType>
double CaMKIIingPosition<BPositionInteractionType>::computeEnergy(double d) {
    
    double U = 0;
    double U_i;
    
    for (auto b: CaMKIIingPoint::getCaMKIIingPoints()) {
        
        Bead* b1 = b->getFirstCylinder()->getFirstBead();
        Bead* b2 = b->getFirstCylinder()->getSecondBead();
        Bead* b3 = b->getSecondCylinder()->getFirstBead();
        
        double kPosition = b->getMCaMKIIingPoint()->getPositionConstant();
        double position = b->getPosition();
        
        if (d == 0.0)
            U_i = _FFType.energy(b1, b2, b3, kPosition, position);
        else
            U_i = _FFType.energy(b1, b2, b3, kPosition, position, d);
        
        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {
            
            //set culprit and return
            _camkiiingCulprit = b;
            
            return -1;
        }
        else
            U += U_i;
    }
    
    return U;
}

template <class BPositionInteractionType>
void CaMKIIingPosition<BPositionInteractionType>::computeForces() {
    
    for (auto b: CaMKIIingPoint::getCaMKIIingPoints()) {
    
        Bead* b1 = b->getFirstCylinder()->getFirstBead();
        Bead* b2 = b->getFirstCylinder()->getSecondBead();
        Bead* b3 = b->getSecondCylinder()->getFirstBead();
        
        double kPosition = b->getMCaMKIIingPoint()->getPositionConstant();
        double position = b->getPosition();
        
        _FFType.forces(b1, b2, b3, kPosition, position);
    }
    
}


template <class BPositionInteractionType>
void CaMKIIingPosition<BPositionInteractionType>::computeForcesAux() {
    
    for (auto b: CaMKIIingPoint::getCaMKIIingPoints()) {
        
        Bead* b1 = b->getFirstCylinder()->getFirstBead();
        Bead* b2 = b->getFirstCylinder()->getSecondBead();
        Bead* b3 = b->getSecondCylinder()->getFirstBead();
        
        double kPosition = b->getMCaMKIIingPoint()->getPositionConstant();
        double position = b->getPosition();
        
        _FFType.forcesAux(b1, b2, b3, kPosition, position);
    }
}


///Template specializations
template double CaMKIIingPosition<CaMKIIingPositionCosine>::computeEnergy(double d);
template void CaMKIIingPosition<CaMKIIingPositionCosine>::computeForces();
template void CaMKIIingPosition<CaMKIIingPositionCosine>::computeForcesAux();