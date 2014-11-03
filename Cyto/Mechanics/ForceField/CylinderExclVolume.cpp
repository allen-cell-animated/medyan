//
//  CylinderExclVolume.cpp
//  Cyto
//
//  Created by Konstantin Popov on 10/31/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "CylinderExclVolume.h"
#include "CylinderExclVolRepulsion.h"
#include "CylinderDB.h"
#include "Cylinder.h"
#include "Bead.h"

template <class CVolumeInteractionType>
double CylinderExclVolume<CVolumeInteractionType>::ComputeEnergy(Cylinder* pc1, Cylinder* pc2, double d) {
    
    Bead* pb1 = pc1->GetFirstBead();
    Bead* pb2 = pc1->GetSecondBead();
    Bead* pb3 = pc1->GetFirstBead();
    Bead* pb4 = pc1->GetSecondBead();
    double kRepuls = 1.0;
    
    
    if (d == 0.0)
        return _FFType.Energy(pb1, pb2, pb3, pb4, kRepuls);
    else
        return _FFType.Energy(pb1, pb2, pb3, pb4, kRepuls, d);   ///This type of function needed for conjugated gradient minimisation only;
    
}

template <class CVolumeInteractionType>
void CylinderExclVolume<CVolumeInteractionType>::ComputeForces(Cylinder* pc1, Cylinder* pc2) {

    Bead* pb1 = pc1->GetFirstBead();
    Bead* pb2 = pc1->GetSecondBead();
    Bead* pb3 = pc1->GetFirstBead();
    Bead* pb4 = pc1->GetSecondBead();
    double kRepuls = 1.0;

    
    _FFType.Forces(pb1, pb2, pb3, pb4, kRepuls );
    
}


template <class CVolumeInteractionType>
void CylinderExclVolume<CVolumeInteractionType>::ComputeForcesAux(Cylinder* pc1, Cylinder* pc2) {
    /// Needed for Conjugated Gradient minimization;
    
    Bead* pb1 = pc1->GetFirstBead();
    Bead* pb2 = pc1->GetSecondBead();
    Bead* pb3 = pc1->GetFirstBead();
    Bead* pb4 = pc1->GetSecondBead();
    double kRepuls = 1.0;
    
    
    _FFType.Forces(pb1, pb2, pb3, pb4, kRepuls );
}

///Template specializations
template double CylinderExclVolume::<CylinderExclVolRepulsion>::ComputeEnergy(Cylinder* pc1, Cylinder* pc2, double d);
template void  CylinderExclVolume::<CylinderExclVolRepulsion>::ComputeForces(Cylinder* pc1, Cylinder* pc2);
template void  CylinderExclVolume::<CylinderExclVolRepulsion>::ComputeForcesAux(Cylinder* pc1, Cylinder* pc2);

