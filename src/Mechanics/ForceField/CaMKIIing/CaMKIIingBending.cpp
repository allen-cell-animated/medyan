
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

#include "CaMKIIingBending.h"

#include "CaMKIIingBendingCosine.h"

#include "CaMKIIingPoint.h"
#include "Cylinder.h"
#include "Bead.h"

template <class BBendingInteractionType>
double CaMKIIingBending<BBendingInteractionType>::computeEnergy(double d) {
    
    double U = 0;
    double U_i;
    
    for (auto b: CaMKIIingPoint::getCaMKIIingPoints()) {
        
        Bead* b1 = b->getFirstCylinder()->getFirstBead();
        Bead* b2 = b->getFirstCylinder()->getSecondBead();
        Bead* b3 = b->getSecondCylinder()->getFirstBead();
        Bead* b4 = b->getSecondCylinder()->getSecondBead();
        double kBend = b->getMCaMKIIingPoint()->getBendingConstant();
        double eqTheta = b->getMCaMKIIingPoint()->getEqTheta();
#if 0
        
        if (d == 0.0)
            U_i = _FFType.energy(b1, b2, b3, b4, kBend, eqTheta);
        else
            U_i = _FFType.energy(b1, b2, b3, b4, kBend, eqTheta, d);
        
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

template <class BBendingInteractionType>
void CaMKIIingBending<BBendingInteractionType>::computeForces() {
   
    for (auto b: CaMKIIingPoint::getCaMKIIingPoints()) {
    
        Bead* b1 = b->getFirstCylinder()->getFirstBead();
        Bead* b2 = b->getFirstCylinder()->getSecondBead();
        Bead* b3 = b->getSecondCylinder()->getFirstBead();
        Bead* b4 = b->getSecondCylinder()->getSecondBead();
        double kBend = b->getMCaMKIIingPoint()->getBendingConstant();
        double eqTheta = b->getMCaMKIIingPoint()->getEqTheta();
        // TODO fix to iterate over each of the _bonds according to coordination number
        //_FFType.forces(b1, b2, b3, b4, kBend, eqTheta);
    }
}

template <class BBendingInteractionType>
void CaMKIIingBending<BBendingInteractionType>::computeForcesAux() {
    
    for (auto b: CaMKIIingPoint::getCaMKIIingPoints()) {
        
        Bead* b1 = b->getFirstCylinder()->getFirstBead();
        Bead* b2 = b->getFirstCylinder()->getSecondBead();
        Bead* b3 = b->getSecondCylinder()->getFirstBead();
        Bead* b4 = b->getSecondCylinder()->getSecondBead();
        double kBend = b->getMCaMKIIingPoint()->getBendingConstant();
        double eqTheta = b->getMCaMKIIingPoint()->getEqTheta();
        // TODO fix to iterate over each of the _bonds according to coordination number
        //_FFType.forcesAux(b1, b2, b3, b4, kBend, eqTheta);
    }
}


///Template specializations
template double CaMKIIingBending<CaMKIIingBendingCosine>::computeEnergy(double d);
template void CaMKIIingBending<CaMKIIingBendingCosine>::computeForces();
template void CaMKIIingBending<CaMKIIingBendingCosine>::computeForcesAux();
