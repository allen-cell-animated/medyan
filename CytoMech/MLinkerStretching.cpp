//
//  MLinkerStratching.cpp
//  Cyto
//
//  Created by Konstantin Popov on 8/28/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "MLinkerStretching.h"

double LinkerStretching::ComputeEnergy(Linker* pl, double d)
{
    if (d == 0.0){
        
        Bead* pb1 = pl->GetFirstCylinder()->GetFirstBead();
        Bead* pb2 = pl->GetFirstCylinder()->GetSecondtBead();
        Bead* pb3 = pl->GetsecondCylinder()->GetFirstBead();
        Bead* pb4 = pl->GetSecondCylinder()->GetSecondtBead();
        double kStretch = pl->GetBendingConst();
        double L = pl->GetEqLength();
        return _FFType.Energy(pb1, pb2, pb3, pb4, kStretch, L);
        
    }
    
    else {
        
        Bead* pb1 = pl->GetFirstCylinder()->GetFirstBead();
        Bead* pb2 = pl->GetFirstCylinder()->GetSecondtBead();
        Bead* pb3 = pl->GetsecondCylinder()->GetFirstBead();
        Bead* pb4 = pl->GetSecondCylinder()->GetSecondtBead();
        double kStretch = pl->GetBendingConst();
        double L = pl->GetEqLength();
        return _FFType.Energy(pb1, pb2, pb3, pb4, kStretch, L, d);   ///This type of function needed for conjugated gradient minimisation only;
    }
}

void LinkerStretching::ComputeForce(Filament* pf)
{
    Bead* pb1 = pl->GetFirstCylinder()->GetFirstBead();
    Bead* pb2 = pl->GetFirstCylinder()->GetSecondtBead();
    Bead* pb3 = pl->GetsecondCylinder()->GetFirstBead();
    Bead* pb4 = pl->GetSecondCylinder()->GetSecondtBead();
    double kStretch = pl->GetBendingConst();
    double L = pl->GetEqLength();

    
        _FFType.Forces(pb1, pb2, pb3, pb4, kStretch, L);
    
}


void LinkerStretching::ComputeForceAux(Filament* pf) /// Needed for Conjugated Gradient minimization;
{
    Bead* pb1 = pl->GetFirstCylinder()->GetFirstBead();
    Bead* pb2 = pl->GetFirstCylinder()->GetSecondtBead();
    Bead* pb3 = pl->GetsecondCylinder()->GetFirstBead();
    Bead* pb4 = pl->GetSecondCylinder()->GetSecondtBead();
    double kStretch = pl->GetBendingConst();
    double L = pl->GetEqLength();

        _FFType.ForcesAux(pb1, pb2, pb3, pb4, kStretch, L);
    
}
