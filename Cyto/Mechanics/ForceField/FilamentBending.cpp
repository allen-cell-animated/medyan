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
double FilamentBending<FBendingInteractionType>::computeEnergy(Filament* f, double d)
{
    if (f->getCylinderVector().size()>1){
    double U = 0.0;
    
        if (d == 0.0){
            for ( auto it = f->getCylinderVector().begin()+1; it != f->getCylinderVector().end(); it++){
                
                auto it2 = it - 1;
                Bead* b1 = (*it2)->getFirstBead();
                Bead* b2 = (*it)->getFirstBead();
                Bead* b3 = (*it)->getSecondBead();
                double k_bend = (*it)->getMCylinder()->getBendingConst();
                
                U += _FFType.energy( b1, b2, b3, k_bend );
            }
        }
        else {
            int index = 0;
            for ( auto it = f->getCylinderVector().begin()+1; it != f->getCylinderVector().end(); it++){
                
                auto it2 = it - 1;
                Bead* b1 = (*it2)->getFirstBead();
                Bead* b2 = (*it)->getFirstBead();
                Bead* b3 = (*it)->getSecondBead();
                double k_bend = (*it)->getMCylinder()->getBendingConst();
                
                U += _FFType.energy( b1, b2, b3, k_bend, d );
                index++;
            }
        }
        //cout << "Bending Energy = " << U << endl;
        
        return U;
    }
    else return 0;
}

template <class FBendingInteractionType>
void FilamentBending<FBendingInteractionType>::computeForces(Filament* f)
{
    
    if (f->getCylinderVector().size()>1){
        for (auto it = f->getCylinderVector().begin()+1; it != f->getCylinderVector().end(); it++){
            
            auto it2 = it - 1;
            Bead* b1 = (*it2)->getFirstBead();
            Bead* b2 = (*it)->getFirstBead();
            Bead* b3 = (*it)->getSecondBead();
            double k_bend = (*it)->getMCylinder()->getBendingConst();
            
            
            _FFType.forces( b1, b2, b3, k_bend );
        }
    }
}

template <class FBendingInteractionType>
void FilamentBending<FBendingInteractionType>::computeForcesAux(Filament* f) /// Needed for Conjugated Gradient minimization;
{
    if (f->getCylinderVector().size()>1){
        for (auto it = f->getCylinderVector().begin()+1; it != f->getCylinderVector().end(); it++){
            
            auto it2 = it - 1;
            Bead* b1 = (*it2)->getFirstBead();
            Bead* b2 = (*it)->getFirstBead();
            Bead* b3 = (*it)->getSecondBead();
            double k_bend = (*it)->getMCylinder()->getBendingConst();

            _FFType.forcesAux( b1, b2, b3, k_bend );
            
        }
    }
}

///Template specializations
template double FilamentBending<FilamentBendingHarmonic>::computeEnergy(Filament* f, double d);
template void  FilamentBending<FilamentBendingHarmonic>::computeForces(Filament* f);
template void  FilamentBending<FilamentBendingHarmonic>::computeForcesAux(Filament* f);
