
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "FilamentBending.h"

#include "FilamentBendingHarmonic.h"
#include "FilamentBendingCosine.h"

#include "Filament.h"
#include "Cylinder.h"

template <class FBendingInteractionType>
double FilamentBending<FBendingInteractionType>::computeEnergy(Filament* f, double d) {
    if (f->getCylinderVector().size()>1){
    double U = 0.0;
    
        if (d == 0.0){
            for ( auto it = f->getCylinderVector().begin()+1;
                       it != f->getCylinderVector().end(); it++){
                
                auto it2 = it - 1;
                Bead* b1 = (*it2)->getFirstBead();
                Bead* b2 = (*it)->getFirstBead();
                Bead* b3 = (*it)->getSecondBead();
                double kBend = (*it)->getMCylinder()->getBendingConst();
                double eqTheta = (*it)->getMCylinder()->getEqTheta();

                U += _FFType.energy(b1, b2, b3, kBend, eqTheta);
            }
        }
        else {
            int index = 0;
            for ( auto it = f->getCylinderVector().begin()+1;
                       it != f->getCylinderVector().end(); it++){
                
                auto it2 = it - 1;
                Bead* b1 = (*it2)->getFirstBead();
                Bead* b2 = (*it)->getFirstBead();
                Bead* b3 = (*it)->getSecondBead();
                double kBend = (*it)->getMCylinder()->getBendingConst();
                double eqTheta = (*it)->getMCylinder()->getEqTheta();
                
                U += _FFType.energy(b1, b2, b3, kBend, eqTheta, d);
                index++;
            }
        }
        return U;
    }
    else return 0;
}

template <class FBendingInteractionType>
void FilamentBending<FBendingInteractionType>::computeForces(Filament* f)
{
    
    if (f->getCylinderVector().size()>1){
        for (auto it = f->getCylinderVector().begin()+1;
                  it != f->getCylinderVector().end(); it++){
            
            auto it2 = it - 1;
            Bead* b1 = (*it2)->getFirstBead();
            Bead* b2 = (*it)->getFirstBead();
            Bead* b3 = (*it)->getSecondBead();
            double kBend = (*it)->getMCylinder()->getBendingConst();
            double eqTheta = (*it)->getMCylinder()->getEqTheta();
            
            _FFType.forces(b1, b2, b3, kBend, eqTheta);
        }
    }
}

template <class FBendingInteractionType>
void FilamentBending<FBendingInteractionType>::computeForcesAux(Filament* f)
{
    if (f->getCylinderVector().size()>1){
        for (auto it = f->getCylinderVector().begin()+1;
                  it != f->getCylinderVector().end(); it++){
            
            auto it2 = it - 1;
            Bead* b1 = (*it2)->getFirstBead();
            Bead* b2 = (*it)->getFirstBead();
            Bead* b3 = (*it)->getSecondBead();
            double kBend = (*it)->getMCylinder()->getBendingConst();
            double eqTheta = (*it)->getMCylinder()->getEqTheta();

            _FFType.forcesAux(b1, b2, b3, kBend, eqTheta);
        }
    }
}

///Template specializations
template double
FilamentBending<FilamentBendingHarmonic>::computeEnergy(Filament* f, double d);
template void
FilamentBending<FilamentBendingHarmonic>::computeForces(Filament* f);
template void
FilamentBending<FilamentBendingHarmonic>::computeForcesAux(Filament* f);
template double
FilamentBending<FilamentBendingCosine>::computeEnergy(Filament* f, double d);
template void
FilamentBending<FilamentBendingCosine>::computeForces(Filament* f);
template void
FilamentBending<FilamentBendingCosine>::computeForcesAux(Filament* f);
