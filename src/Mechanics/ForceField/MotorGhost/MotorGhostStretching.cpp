
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

#include "MotorGhostStretching.h"

#include "MotorGhostStretchingHarmonic.h"

#include "MotorGhost.h"
#include "Cylinder.h"
#include "Bead.h"
#include "cross_check.h"

template <class MStretchingInteractionType>
void MotorGhostStretching<MStretchingInteractionType>::vectorize() {
    
    beadSet = new int[n * MotorGhost::getMotorGhosts().size()];
    kstr = new double[MotorGhost::getMotorGhosts().size()];
    eql = new double[MotorGhost::getMotorGhosts().size()];
    pos1 = new double[MotorGhost::getMotorGhosts().size()];
    pos2 = new double[MotorGhost::getMotorGhosts().size()];
    
    int i = 0;
    
    for (auto m: MotorGhost::getMotorGhosts()) {
        
        beadSet[n * i] = m->getFirstCylinder()->getFirstBead()->_dbIndex;
        beadSet[n * i + 1] = m->getFirstCylinder()->getSecondBead()->_dbIndex;
        beadSet[n * i + 2] = m->getSecondCylinder()->getFirstBead()->_dbIndex;
        beadSet[n * i + 3] = m->getSecondCylinder()->getSecondBead()->_dbIndex;
        
        kstr[i] = m->getMMotorGhost()->getStretchingConstant();
        eql[i] = m->getMMotorGhost()->getEqLength();
        pos1[i] = m->getFirstPosition();
        pos2[i] = m->getSecondPosition();
        
        i++;
    }
}

template<class MStretchingInteractionType>
void MotorGhostStretching<MStretchingInteractionType>::deallocate() {
    
    delete beadSet;
    delete kstr;
    delete eql;
    delete pos1;
    delete pos2;
}


template <class MStretchingInteractionType>
double MotorGhostStretching<MStretchingInteractionType>::computeEnergy(double* coord, double *f, double d){
    
    double U_i;
    
    if (d == 0.0)
        U_i = _FFType.energy(coord, f, beadSet, kstr, eql, pos1, pos2);
    else
        U_i = _FFType.energy(coord, f, beadSet, kstr, eql, pos1, pos2, d);
    
#ifdef CROSSCHECK
    double U2 = 0;
    double U_ii;
    
    for (auto m: MotorGhost::getMotorGhosts()) {
        
        Bead* b1 = m->getFirstCylinder()->getFirstBead();
        Bead* b2 = m->getFirstCylinder()->getSecondBead();
        Bead* b3 = m->getSecondCylinder()->getFirstBead();
        Bead* b4 = m->getSecondCylinder()->getSecondBead();
        
        double kStretch = m->getMMotorGhost()->getStretchingConstant();
        double eqLength = m->getMMotorGhost()->getEqLength();
        
        double pos1 = m->getFirstPosition();
        double pos2 = m->getSecondPosition();
        
        if (d == 0.0)
            U_ii = _FFType.energy(b1, b2, b3, b4, pos1, pos2, kStretch, eqLength);
        else
            U_ii = _FFType.energy(b1, b2, b3, b4, pos1, pos2, kStretch, eqLength, d);
        
        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_ii != U_ii || U_ii < -1.0) {
            
            //set culprit and return
            _motorCulprit = m;
            
            U2=-1;
            break;
        }
        else
            U2 += U_ii;
    }
    if(abs(U_i-U2)<=U2/100000000000)
        std::cout<<"E M YES "<<endl;
    else
    {   std::cout<<U_i<<" "<<U2<<endl;
        exit(EXIT_FAILURE);
    }

    
#endif
    
    return U_i;
}

template <class MStretchingInteractionType>
void MotorGhostStretching<MStretchingInteractionType>::computeForces(double *coord, double *f) {
    
    _FFType.forces(coord, f, beadSet, kstr, eql, pos1, pos2);
#ifdef CROSSCHECK
    for (auto m: MotorGhost::getMotorGhosts()) {
        
        Bead* b1 = m->getFirstCylinder()->getFirstBead();
        Bead* b2 = m->getFirstCylinder()->getSecondBead();
        Bead* b3 = m->getSecondCylinder()->getFirstBead();
        Bead* b4 = m->getSecondCylinder()->getSecondBead();
        double kStretch = m->getMMotorGhost()->getStretchingConstant();
        double eqLength = m->getMMotorGhost()->getEqLength();
        
        double pos1 = m->getFirstPosition();
        double pos2 = m->getSecondPosition();
        
        if(cross_checkclass::Aux)
        {double f0 = _FFType.forcesAux(b1, b2, b3, b4, pos1, pos2, kStretch, eqLength);
//            m->getMMotorGhost()->stretchForce = f0;
        }
        else
        {double f0 = _FFType.forces(b1, b2, b3, b4, pos1, pos2, kStretch, eqLength);
//            m->getMMotorGhost()->stretchForce = f0;
        }
    }
    if(cross_checkclass::Aux){
        auto state=cross_check::crosscheckAuxforces(f);
        std::cout<<"F S+B+L+M YES "<<state<<endl;}
    else{
        auto state=cross_check::crosscheckforces(f);
        std::cout<<"F S+B+L+M YES "<<state<<endl;}
#endif
}


///Temlate specializations
template double MotorGhostStretching<MotorGhostStretchingHarmonic>::computeEnergy(double *coord, double *f, double d);
template void MotorGhostStretching<MotorGhostStretchingHarmonic>::computeForces(double *coord, double *f);
template void MotorGhostStretching<MotorGhostStretchingHarmonic>::vectorize();
template void MotorGhostStretching<MotorGhostStretchingHarmonic>::deallocate();


