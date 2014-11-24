//
//  CylinderExclVolume.cpp
//  Cyto
//
//  Created by Konstantin Popov on 10/31/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "CylinderExclVolume.h"

#include "CylinderExclVolRepulsion.h"

#include "Cylinder.h"
#include "Bead.h"

template <class CVolumeInteractionType>
double CylinderExclVolume<CVolumeInteractionType>::computeEnergy(Cylinder* c1, Cylinder* c2, double d) {
    
    Bead* b1 = c1->getFirstBead();
    Bead* b2 = c1->getSecondBead();
    Bead* b3 = c2->getFirstBead();
    Bead* b4 = c2->getSecondBead();
    double kRepuls = c1->getMCylinder()->getExVolConst();
    
    
    if (d == 0.0)
        return _FFType.energy(b1, b2, b3, b4, kRepuls);
    else
        return _FFType.energy(b1, b2, b3, b4, kRepuls, d);   ///This type of function needed for conjugated gradient minimisation only;
    
}

template <class CVolumeInteractionType>
void CylinderExclVolume<CVolumeInteractionType>::computeForces(Cylinder* c1, Cylinder* c2) {

    Bead* b1 = c1->getFirstBead();
    Bead* b2 = c1->getSecondBead();
    Bead* b3 = c2->getFirstBead();
    Bead* b4 = c2->getSecondBead();
    double kRepuls = c1->getMCylinder()->getExVolConst();

    
    _FFType.forces(b1, b2, b3, b4, kRepuls );
    
}


template <class CVolumeInteractionType>
void CylinderExclVolume<CVolumeInteractionType>::computeForcesAux(Cylinder* c1, Cylinder* c2) {
    /// Needed for Conjugated Gradient minimization;
    
    Bead* b1 = c1->getFirstBead();
    Bead* b2 = c1->getSecondBead();
    Bead* b3 = c2->getFirstBead();
    Bead* b4 = c2->getSecondBead();
    double kRepuls = c1->getMCylinder()->getExVolConst();
    
    
    _FFType.forcesAux(b1, b2, b3, b4, kRepuls );
}

///Template specializations
template double CylinderExclVolume<CylinderExclVolRepulsion>::computeEnergy(Cylinder* c1, Cylinder* c2, double d);
template void  CylinderExclVolume<CylinderExclVolRepulsion>::computeForces(Cylinder* c1, Cylinder* c2);
template void  CylinderExclVolume<CylinderExclVolRepulsion>::computeForcesAux(Cylinder* c1, Cylinder* c2);

