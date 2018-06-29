
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

#include "MotorGhostStretching.h"

#include "MotorGhostStretchingHarmonic.h"

#include "MotorGhost.h"
#include "Cylinder.h"
#include "Bead.h"

template <class MStretchingInteractionType>
double MotorGhostStretching<MStretchingInteractionType>::computeEnergy(double d) {
    
    double U = 0;
    double U_i;
    
    for (auto m: MotorGhost::getMotorGhosts()) {
        
        Bead* b1 = m->getFirstCylinder()->getFirstBead();
        Bead* b2 = m->getFirstCylinder()->getSecondBead();
        Bead* b3 = m->getSecondCylinder()->getFirstBead();
        Bead* b4 = m->getSecondCylinder()->getSecondBead();
        double kStretch = m->getMMotorGhost()->getStretchingConstant();
        double eqLength = m->getMMotorGhost()->getEqLength();
        
        double pos1 = m->getFirstPosition();
        double pos2 = m->getSecondPosition();
        
        if (d == 0.0)
            U_i = _FFType.energy(b1, b2, b3, b4, pos1, pos2, kStretch, eqLength);
        else
            U_i = _FFType.energy(b1, b2, b3, b4, pos1, pos2, kStretch, eqLength, d);
        
        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {
            
            //set culprit and return
            _motorCulprit = m;
            
            return -1;
        }
        else
            U += U_i;
    }
    return U;
}

template <class MStretchingInteractionType>
void MotorGhostStretching<MStretchingInteractionType>::computeForces() {
    
    for (auto m: MotorGhost::getMotorGhosts()) {
    
        Bead* b1 = m->getFirstCylinder()->getFirstBead();
        Bead* b2 = m->getFirstCylinder()->getSecondBead();
        Bead* b3 = m->getSecondCylinder()->getFirstBead();
        Bead* b4 = m->getSecondCylinder()->getSecondBead();
        double kStretch = m->getMMotorGhost()->getStretchingConstant();
        double eqLength = m->getMMotorGhost()->getEqLength();
        
        double pos1 = m->getFirstPosition();
        double pos2 = m->getSecondPosition();
        
        double f0 = _FFType.forces(b1, b2, b3, b4, pos1, pos2, kStretch, eqLength);
        m->getMMotorGhost()->stretchForce = f0;
    }
    double maxF = 0;
    
    //calc max force
    for(auto b: Bead::getBeads())
        maxF = max(maxF, sqrt(b->FADotFA()));
    std::cout<<"maxF "<<getName()<<" "<<maxF<<endl;
}


template <class MStretchingInteractionType>
void MotorGhostStretching<MStretchingInteractionType>::computeForcesAux() {

    for (auto m: MotorGhost::getMotorGhosts()) {
    
        Bead* b1 = m->getFirstCylinder()->getFirstBead();
        Bead* b2 = m->getFirstCylinder()->getSecondBead();
        Bead* b3 = m->getSecondCylinder()->getFirstBead();
        Bead* b4 = m->getSecondCylinder()->getSecondBead();
        double kStretch = m->getMMotorGhost()->getStretchingConstant();
        double eqLength = m->getMMotorGhost()->getEqLength();
        
        double pos1 = m->getFirstPosition();
        double pos2 = m->getSecondPosition();
        
        double f0 = _FFType.forcesAux(b1, b2, b3, b4, pos1, pos2, kStretch, eqLength);
        m->getMMotorGhost()->stretchForce = f0;
    }
    double maxF = 0;
    
    //calc max force
    for(auto b: Bead::getBeads())
        maxF = max(maxF, sqrt(b->FADotFA()));
    std::cout<<"maxF "<<getName()<<" "<<maxF<<endl;
}

///Template specializations
template double MotorGhostStretching<MotorGhostStretchingHarmonic>::computeEnergy(double d);
template void MotorGhostStretching<MotorGhostStretchingHarmonic>::computeForces();
template void MotorGhostStretching<MotorGhostStretchingHarmonic>::computeForcesAux();

