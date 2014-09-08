//
//  MMotorGhostStretching.cpp
//  Cyto
//
//  Created by Konstantin Popov on 9/3/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "MMotorGhostStretching.h"
template <class MStretchingInteractionType>
double MotorGhostStretching<MStretchingInteractionType>::ComputeEnergy(MotorGhost* pm, double d)
{
    if (d == 0.0){
        
        Bead* pb1 = pm->GetFirstCylinder()->GetFirstBead();
        Bead* pb2 = pm->GetFirstCylinder()->GetSecondtBead();
        Bead* pb3 = pm->GetsecondCylinder()->GetFirstBead();
        Bead* pb4 = pm->GetSecondCylinder()->GetSecondtBead();
        double kStretch = pm->GetBendingConst();
        double L = pm->GetEqLength();
        return _FFType.Energy(pb1, pb2, pb3, pb4, kStretch, L);
        
    }
    
    else {
        
        Bead* pb1 = pm->GetFirstCylinder()->GetFirstBead();
        Bead* pb2 = pm->GetFirstCylinder()->GetSecondtBead();
        Bead* pb3 = pm->GetsecondCylinder()->GetFirstBead();
        Bead* pb4 = pm->GetSecondCylinder()->GetSecondtBead();
        double kStretch = pm->GetBendingConst();
        double L = pm->GetEqLength();

        return _FFType.Energy(pb1, pb2, pb3, pb4, kStretch, L, d);   ///This type of function needed for conjugated gradient minimisation only;
    }
}

template <class MStretchingInteractionType>
void MotorGhostStretching<MStretchingInteractionType>::ComputeForces(MotorGhost* pm)
{
    Bead* pb1 = pm->GetFirstCylinder()->GetFirstBead();
    Bead* pb2 = pm->GetFirstCylinder()->GetSecondtBead();
    Bead* pb3 = pm->GetsecondCylinder()->GetFirstBead();
    Bead* pb4 = pm->GetSecondCylinder()->GetSecondtBead();
    double kStretch = pm->GetBendingConst();
    double L = pm->GetEqLength();

    
    
    _FFType.Forces(pb1, pb2, pb3, pb4, kStretch, L);
    
}


template <class MStretchingInteractionType>
void MotorGhostStretching<MStretchingInteractionType>::ComputeForcesAux(MotorGhost* pm) /// Needed for Conjugated Gradient minimization;
{
    Bead* pb1 = pm->GetFirstCylinder()->GetFirstBead();
    Bead* pb2 = pm->GetFirstCylinder()->GetSecondtBead();
    Bead* pb3 = pm->GetsecondCylinder()->GetFirstBead();
    Bead* pb4 = pm->GetSecondCylinder()->GetSecondtBead();
    double kStretch = pm->GetBendingConst();
    double L = pm->GetEqLength();
    
    _FFType.ForcesAux(pb1, pb2, pb3, pb4, kStretch, L);
    
}
