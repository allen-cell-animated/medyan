#ifdef CAMKII
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "CamkiiBending.h"

#include "CamkiiBendingCosine.h"

#include "Camkii.h"
#include "Cylinder.h"

template <class FBendingInteractionType>
double CamkiiBending<FBendingInteractionType>::computeEnergy(double d) {
    
    double U = 0;
    double U_i;
    
    for (auto f: Camkii::getCamkiis()) {
        
        U_i = 0;
        
        if (f->getCylinders().size() > 1){
            
            if (d == 0.0){
                for (auto it = f->getCylinders().begin()+1;
                          it != f->getCylinders().end(); it++){
                    
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
                for (auto it = f->getCylinders().begin()+1;
                          it != f->getCylinders().end(); it++){
                    
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
            _camkiiCulprit = f;
            
            return -1;
        }
        else
            U += U_i;
    }
    
    return U;
}

template <class FBendingInteractionType>
void CamkiiBending<FBendingInteractionType>::computeForces()
{
    for (auto f: Camkii::getCamkiis()) {
        
        if (f->getCylinders().size()>1){
            for (auto it = f->getCylinders().begin()+1;
                      it != f->getCylinders().end(); it++){
                
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

template <class FBendingInteractionType>
void CamkiiBending<FBendingInteractionType>::computeForcesAux()
{
    for (auto f: Camkii::getCamkiis()) {
        
        if (f->getCylinders().size()>1){
            for (auto it = f->getCylinders().begin()+1;
                 it != f->getCylinders().end(); it++){
                
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
//template double CamkiiBending<CamkiiBendingHarmonic>::computeEnergy(double d);
//template void CamkiiBending<CamkiiBendingHarmonic>::computeForces();
//template void CamkiiBending<CamkiiBendingHarmonic>::computeForcesAux();
template double CamkiiBending<CamkiiBendingCosine>::computeEnergy(double d);
template void CamkiiBending<CamkiiBendingCosine>::computeForces();
template void CamkiiBending<CamkiiBendingCosine>::computeForcesAux();
#endif