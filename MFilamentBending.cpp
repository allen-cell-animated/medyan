//
//  MFilamentBending.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "MFilamentBending.h"

double FilamentBending::ComputeEnergy(Filament* pf, double d)
{
    double U = 0.0;
    
    if (d == 0.0){
        for ( auto it = pf->_pCylinderVector.begin()+1; it != _pCylinderVector.end(); it++){
            
            Bead* pb1 = (it-1)->GetFirstBead();
            Bead* pb2 = it->GetFirstBead();
            Bead* pb3 = it->GetSecondBead();
            double k_bend = it->GetBendingConst();
            
            U += _FFType.Energy( pb1, pb2, pb3, k_bend );
            return U;
        }
    }
    else {
        for ( auto it = pf->_pCylinderVector.begin()+1; it != _pCylinderVector.end(); it++){
            
            Bead* pb1 = (it-1)->GetFirstBead();
            Bead* pb2 = it->GetFirstBead();
            Bead* pb3 = it->GetSecondBead();
            double k_bend = it->GetBendingConst();
            
            U += _FFType.Energy( pb1, pb2, pb3, k_bend, d );
            return U;
        }    }
}

void FilamentBending::ComputeForce(Filament* pf)
{
    for ( auto it = _pCylinderVector.begin()+1; it != _pCylinderVector.end(); it++){
        
        Bead* pb1 = (it-1)->GetFirstBead();
        Bead* pb2 = it->GetFirstBead();
        Bead* pb3 = it->GetSecondBead();
        double k_bend = it->GetBendingConst();
        
        _FFType.Forces( pb1, pb2, pb3, k_bend );
    
    }
    
}


void FilamentBending::ComputeForceAux(Filament* pf) /// Needed for Conjugated Gradient minimization;
{
    for ( auto it = _pCylinderVector.begin()+1; it != _pCylinderVector.end(); it++){
        
        
        Bead* pb1 = (it-1)->GetFirstBead();
        Bead* pb2 = it->GetFirstBead();
        Bead* pb3 = it->GetSecondBead();
        double k_bend = it->GetBendingConst();
        
        _FFType.ForcesAux( pb1, pb2, pb3, k_bend );
        
    }
}