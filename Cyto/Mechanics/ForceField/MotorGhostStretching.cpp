//
//  MotorGhostStretching.cpp
//  Cyto
//
//  Created by Konstantin Popov on 9/3/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "MotorGhostStretching.h"
#include "MotorGhostStretchingHarmonic.h"
#include "MotorGhostDB.h"
#include "Cylinder.h"
#include "Bead.h"

template <class MStretchingInteractionType>
double MotorGhostStretching<MStretchingInteractionType>::computeEnergy(MotorGhost* pm, double d) {
    
    Bead* b1 = pm->getFirstCylinder()->getFirstBead();
    Bead* b2 = pm->getFirstCylinder()->getSecondBead();
    Bead* b3 = pm->getSecondCylinder()->getFirstBead();
    Bead* b4 = pm->getSecondCylinder()->getSecondBead();
    double kStretch = pm->getMMotorGhost()->getStretchingConstant();
    double L = pm->getMMotorGhost()->getEqLength();
    
    double pos1 = pm->getFirstPosition();
    double pos2 = pm->getSecondPosition();
    
    if (d == 0.0)
        return _FFType.energy(b1, b2, b3, b4, pos1, pos2, kStretch, L);
    else
        return _FFType.energy(b1, b2, b3, b4, pos1, pos2, kStretch, L, d);   ///This type of function needed for conjugated gradient minimisation only;
    
}

template <class MStretchingInteractionType>
void MotorGhostStretching<MStretchingInteractionType>::computeForces(MotorGhost* pm) {
    
    Bead* b1 = pm->getFirstCylinder()->getFirstBead();
    Bead* b2 = pm->getFirstCylinder()->getSecondBead();
    Bead* b3 = pm->getSecondCylinder()->getFirstBead();
    Bead* b4 = pm->getSecondCylinder()->getSecondBead();
    double kStretch = pm->getMMotorGhost()->getStretchingConstant();
    double L = pm->getMMotorGhost()->getEqLength();
    
    double pos1 = pm->getFirstPosition();
    double pos2 = pm->getSecondPosition();
    
    _FFType.forces(b1, b2, b3, b4, pos1, pos2, kStretch, L);
    
}


template <class MStretchingInteractionType>
void MotorGhostStretching<MStretchingInteractionType>::computeForcesAux(MotorGhost* pm) {/// Needed for Conjugated Gradient minimization;

    Bead* b1 = pm->getFirstCylinder()->getFirstBead();
    Bead* b2 = pm->getFirstCylinder()->getSecondBead();
    Bead* b3 = pm->getSecondCylinder()->getFirstBead();
    Bead* b4 = pm->getSecondCylinder()->getSecondBead();
    double kStretch = pm->getMMotorGhost()->getStretchingConstant();
    double L = pm->getMMotorGhost()->getEqLength();
    
    double pos1 = pm->getFirstPosition();
    double pos2 = pm->getSecondPosition();
    
    _FFType.forcesAux(b1, b2, b3, b4, pos1, pos2, kStretch, L);
}

///Template specializations
template double MotorGhostStretching<MotorGhostStretchingHarmonic>::computeEnergy(MotorGhost* pm, double d);
template void  MotorGhostStretching<MotorGhostStretchingHarmonic>::computeForces(MotorGhost* pm);
template void  MotorGhostStretching<MotorGhostStretchingHarmonic>::computeForcesAux(MotorGhost* pm);

