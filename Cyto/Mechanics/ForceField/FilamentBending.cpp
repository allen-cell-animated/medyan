//
//  FilamentBending.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "FilamentBending.h"
#include "FilamentBendingHarmonic.h"
#include "Filament.h"

template <class FBendingInteractionType>
double FilamentBending<FBendingInteractionType>::ComputeEnergy(Filament* pf, double d)
{
    if (pf->getCylinderVector().size()>1){
    double U = 0.0;
    
        if (d == 0.0){
            for ( auto it = pf->getCylinderVector().begin()+1; it != pf->getCylinderVector().end(); it++){
                
                auto it2 = it - 1;
                Bead* pb1 = (*it2)->GetFirstBead();
                Bead* pb2 = (*it)->GetFirstBead();
                Bead* pb3 = (*it)->GetSecondBead();
                double k_bend = (*it)->getMCylinder()->GetBendingConst();
                
                U += _FFType.Energy( pb1, pb2, pb3, k_bend );
            }
        }
        else {
            int index = 0;
            for ( auto it = pf->getCylinderVector().begin()+1; it != pf->getCylinderVector().end(); it++){
                
                auto it2 = it - 1;
                Bead* pb1 = (*it2)->GetFirstBead();
                Bead* pb2 = (*it)->GetFirstBead();
                Bead* pb3 = (*it)->GetSecondBead();
                double k_bend = (*it)->getMCylinder()->GetBendingConst();
                
                U += _FFType.Energy( pb1, pb2, pb3, k_bend, d );
                index++;
            }
        }
        //std::cout << "Bending Energy = " << U << std::endl;
        
        return U;
    }
    else return 0;
}

template <class FBendingInteractionType>
void FilamentBending<FBendingInteractionType>::ComputeForces(Filament* pf)
{
    
    if (pf->getCylinderVector().size()>1){
        for (auto it = pf->getCylinderVector().begin()+1; it != pf->getCylinderVector().end(); it++){
            
            auto it2 = it - 1;
            Bead* pb1 = (*it2)->GetFirstBead();
            Bead* pb2 = (*it)->GetFirstBead();
            Bead* pb3 = (*it)->GetSecondBead();
            double k_bend = (*it)->getMCylinder()->GetBendingConst();
            
            
            _FFType.Forces( pb1, pb2, pb3, k_bend );
        }
    }
}

template <class FBendingInteractionType>
void FilamentBending<FBendingInteractionType>::ComputeForcesAux(Filament* pf) /// Needed for Conjugated Gradient minimization;
{
    if (pf->getCylinderVector().size()>1){
        for (auto it = pf->getCylinderVector().begin()+1; it != pf->getCylinderVector().end(); it++){
            
            auto it2 = it - 1;
            Bead* pb1 = (*it2)->GetFirstBead();
            Bead* pb2 = (*it)->GetFirstBead();
            Bead* pb3 = (*it)->GetSecondBead();
            double k_bend = (*it)->getMCylinder()->GetBendingConst();

            _FFType.ForcesAux( pb1, pb2, pb3, k_bend );
            
        }
    }
}

///Template specializations
template double FilamentBending<FilamentBendingHarmonic>::ComputeEnergy(Filament* pf, double d);
template void  FilamentBending<FilamentBendingHarmonic>::ComputeForces(Filament* pf);
template void  FilamentBending<FilamentBendingHarmonic>::ComputeForcesAux(Filament* pf);
