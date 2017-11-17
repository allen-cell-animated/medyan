/*

Implementing the "easy" bending force field:

TODO: Implement this

*/

#include "MembraneBending.h"

// TODO: change this
#include "FilamentBendingHarmonic.h"
#include "FilamentBendingCosine.h"

#include "Membrane.h"
#include "Edge.h"

template <class MembraneBendingInteractionType>
double MembraneBending<MembraneBendingInteractionType>::computeEnergy(double d) {
    
    double U = 0;
    double U_i;

    // TODO: all the following
    
    for (auto f: Filament::getFilaments()) {
        
        U_i = 0;
        
        if (f->getCylinderVector().size() > 1){
            
            if (d == 0.0){
                for (auto it = f->getCylinderVector().begin()+1;
                          it != f->getCylinderVector().end(); it++){
                    
                    auto it2 = it - 1;
                    Bead* b1 = (*it2)->getFirstBead();
                    Bead* b2 = (*it)->getFirstBead();
                    Bead* b3 = (*it)->getSecondBead();
                    double kBend = (*it)->getMCylinder()->getBendingConst();
                    double eqTheta = (*it)->getMCylinder()->getEqTheta();
                    
                    U_i += _FFType.energy(b1, b2, b3, kBend, eqTheta);
                }
            }
            else {
                for (auto it = f->getCylinderVector().begin()+1;
                          it != f->getCylinderVector().end(); it++){
                    
                    auto it2 = it - 1;
                    Bead* b1 = (*it2)->getFirstBead();
                    Bead* b2 = (*it)->getFirstBead();
                    Bead* b3 = (*it)->getSecondBead();
                    double kBend = (*it)->getMCylinder()->getBendingConst();
                    double eqTheta = (*it)->getMCylinder()->getEqTheta();
                    
                    U_i += _FFType.energy(b1, b2, b3, kBend, eqTheta, d);
                }
            }
        }
        
        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {
            
            //set culprit and return
            _filamentCulprit = f;
            
            return -1;
        }
        else
            U += U_i;
    }
    
    return U;
}

template <class MembraneBendingInteractionType>
void MembraneBending<MembraneBendingInteractionType>::computeForces()
{
    for (auto f: Filament::getFilaments()) {
        
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
}

template <class MembraneBendingInteractionType>
void MembraneBending<MembraneBendingInteractionType>::computeForcesAux()
{
    for (auto f: Filament::getFilaments()) {
        
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
}

///Template specializations
template double FilamentBending<FilamentBendingHarmonic>::computeEnergy(double d);
template void FilamentBending<FilamentBendingHarmonic>::computeForces();
template void FilamentBending<FilamentBendingHarmonic>::computeForcesAux();
template double FilamentBending<FilamentBendingCosine>::computeEnergy(double d);
template void FilamentBending<FilamentBendingCosine>::computeForces();
template void FilamentBending<FilamentBendingCosine>::computeForcesAux();
