//
//  LinkerStratching.cpp
//  Cyto
//
//  Created by Konstantin Popov on 8/28/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "LinkerStretching.h"
#include "LinkerStretchingHarmonic.h"
#include "Cylinder.h"
#include "Bead.h"
#include "Linker.h"

template <class LStretchingInteractionType>
double LinkerStretching<LStretchingInteractionType>::ComputeEnergy(Linker* pl, double d) {
    Bead* pb1 = pl->GetFirstCylinder()->getMCylinder()->GetFirstBead();
    Bead* pb2 = pl->GetFirstCylinder()->getMCylinder()->GetSecondBead();
    Bead* pb3 = pl->GetSecondCylinder()->getMCylinder()->GetFirstBead();
    Bead* pb4 = pl->GetSecondCylinder()->getMCylinder()->GetSecondBead();
    double kStretch = pl->GetStretchingConstant();
    double L = pl->GetEqLength();
    double pos1 = pl->GetFirstPosition();
    double pos2 = pl->GetSecondPosition();
    
    
    if (d == 0.0)
        return _FFType.Energy(pb1, pb2, pb3, pb4, pos1, pos2, kStretch, L);
    else
        return _FFType.Energy(pb1, pb2, pb3, pb4, pos1, pos2, kStretch, L, d);   ///This type of function needed for conjugated gradient minimisation only;
    
}

template <class LStretchingInteractionType>
void LinkerStretching<LStretchingInteractionType>::ComputeForces(Linker* pl) {
    Bead* pb1 = pl->GetFirstCylinder()->getMCylinder()->GetFirstBead();
    Bead* pb2 = pl->GetFirstCylinder()->getMCylinder()->GetSecondBead();
    Bead* pb3 = pl->GetSecondCylinder()->getMCylinder()->GetFirstBead();
    Bead* pb4 = pl->GetSecondCylinder()->getMCylinder()->GetSecondBead();
    double kStretch = pl->GetStretchingConstant();
    double L = pl->GetEqLength();
    
    double pos1 = pl->GetFirstPosition();
    double pos2 = pl->GetSecondPosition();
    
    _FFType.Forces(pb1, pb2, pb3, pb4, pos1, pos2, kStretch, L);
    
}


template <class LStretchingInteractionType>
void LinkerStretching<LStretchingInteractionType>::ComputeForcesAux(Linker* pl) { /// Needed for Conjugated Gradient minimization;

    Bead* pb1 = pl->GetFirstCylinder()->getMCylinder()->GetFirstBead();
    Bead* pb2 = pl->GetFirstCylinder()->getMCylinder()->GetSecondBead();
    Bead* pb3 = pl->GetSecondCylinder()->getMCylinder()->GetFirstBead();
    Bead* pb4 = pl->GetSecondCylinder()->getMCylinder()->GetSecondBead();
    double kStretch = pl->GetStretchingConstant();
    double L = pl->GetEqLength();
    
    double pos1 = pl->GetFirstPosition();
    double pos2 = pl->GetSecondPosition();
    
    _FFType.ForcesAux(pb1, pb2, pb3, pb4, pos1, pos2, kStretch, L);
    
}

///Template specializations
template double LinkerStretching<LinkerStretchingHarmonic>::ComputeEnergy(Linker* pl, double d);
template void  LinkerStretching<LinkerStretchingHarmonic>::ComputeForces(Linker* pl);
template void  LinkerStretching<LinkerStretchingHarmonic>::ComputeForcesAux(Linker* pl);
