//
//  FilamentBending.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "FilamentBending.h"

template <class FBendingInteractionType>
double FilamentBending<FBendingInteractionType>::ComputeEnergy(Filament* pf, double d)
{
    double U = 0.0;
    
    if (d == 0.0){
        for ( auto it = pf->getCylinderVector().begin()+1; it != pf->getCylinderVector().end(); it++){
            
            Bead* pb1 = (*(it-1))->getMCylinder()->GetFirstBead();
            Bead* pb2 = (*it)->getMCylinder()->GetFirstBead();
            Bead* pb3 = (*it)->getMCylinder()->GetSecondBead();
            double k_bend = (*it)->getMCylinder()->GetBendingConst();
            
            U += _FFType.Energy( pb1, pb2, pb3, k_bend );
            return U;
        }
    }
    else {
        for ( auto it = pf->getCylinderVector().begin()+1; it != pf->getCylinderVector().end(); it++){
            
            Bead* pb1 = (*(it-1))->getMCylinder()->GetFirstBead();
            Bead* pb2 = (*it)->getMCylinder()->GetFirstBead();
            Bead* pb3 = (*it)->getMCylinder()->GetSecondBead();
            double k_bend = (*it)->getMCylinder()->GetBendingConst();
            
            U += _FFType.Energy( pb1, pb2, pb3, k_bend, d );
            return U;
        }    }
}
template <class FBendingInteractionType>
void FilamentBending<FBendingInteractionType>::ComputeForces(Filament* pf)
{
    for (auto it = pf->getCylinderVector().begin()+1; it != pf->getCylinderVector().end(); it++){
        
        Bead* pb1 = (*(it-1))->getMCylinder()->GetFirstBead();
        Bead* pb2 = (*it)->getMCylinder()->GetFirstBead();
        Bead* pb3 = (*it)->getMCylinder()->GetSecondBead();
        double k_bend = (*it)->getMCylinder()->GetBendingConst();
        
        
        _FFType.Forces( pb1, pb2, pb3, k_bend );
    
    }
    
}

template <class FBendingInteractionType>
void FilamentBending<FBendingInteractionType>::ComputeForcesAux(Filament* pf) /// Needed for Conjugated Gradient minimization;
{
    for (auto it = pf->getCylinderVector().begin()+1; it != pf->getCylinderVector().end(); it++){
        
        Bead* pb1 = (*(it-1))->getMCylinder()->GetFirstBead();
        Bead* pb2 = (*it)->getMCylinder()->GetFirstBead();
        Bead* pb3 = (*it)->getMCylinder()->GetSecondBead();
        double k_bend = (*it)->getMCylinder()->GetBendingConst();
        
        
        _FFType.ForcesAux( pb1, pb2, pb3, k_bend );
        
    }
}