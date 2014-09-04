//
//  MFilamentStretching.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "MFilamentStratching.h"

double FilamentStretching::ComputeEnergy(Filament* pf, double d)
{
    double U = 0.0;
    
    if (d == 0.0){
        for(auto it : pf->getCylinderVector()){
            
            Bead* pb1 = it->GetFirstBead();
            Bead* pb2 = it->GetFirstBead();
            double kStr = it->GetStretchingConst();
            double L = it->GetEqLength();
            U += _FFType.Energy(pb1, pb2, kStr, L);
            return U;
        }
    }
    else {
        for(auto it : pf->getCylinderVector()){
            Bead* pb1 = it->GetFirstBead();
            Bead* pb2 = it->GetFirstBead();
            double kStr =it->GetStretchingConst();
            double L = it->GetEqLength();
            
            U += _FFType.Energy(pb1, pb2, k_str, L, d);   ///This type of function needed for conjugated gradient minimisation only;
            return U;
        }
    }
}

void FilamentStretching::ComputeForce(Filament* pf)
{
   for(auto it : pf->getCylinderVector()){
       
       Bead* pb1 = it->GetFirstBead();
       Bead* pb2 = it->GetFirstBead();
       double kStr =it->GetStretchingConst();
       double L = it->GetEqLength();
       _FFType.Forces(pb1, pb2, k_str, L);
   }
}


void FilamentStretching::ComputeForceAux(Filament* pf) /// Needed for Conjugated Gradient minimization;
{
    for(auto it : pf->getCylinderVector()){
        
        Bead* pb1 = it->GetFirstBead();
        Bead* pb2 = it->GetFirstBead();
        double kStr =it->GetStretchingConst();
        double L = it->GetEqLength();
        
        _FFType.ForcesAux(pb1, pb2, kStr, L);
    }
}
