
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

#include "LinkerStretching.h"

#include "LinkerStretchingHarmonic.h"

#include "Cylinder.h"
#include "Linker.h"
#include "Bead.h"
#include "cross_check.h"
template <class LStretchingInteractionType>
void LinkerStretching<LStretchingInteractionType>::vectorize() {
    
    beadSet = new int[n * Linker::getLinkers().size()];
    kstr = new double[Linker::getLinkers().size()];
    eql = new double[Linker::getLinkers().size()];
    pos1 = new double[Linker::getLinkers().size()];
    pos2 = new double[Linker::getLinkers().size()];
    
    int i = 0;
    
    for (auto l: Linker::getLinkers()) {
        
        beadSet[n * i] = l->getFirstCylinder()->getFirstBead()->_dbIndex;
        beadSet[n * i + 1] = l->getFirstCylinder()->getSecondBead()->_dbIndex;
        beadSet[n * i + 2] = l->getSecondCylinder()->getFirstBead()->_dbIndex;
        beadSet[n * i + 3] = l->getSecondCylinder()->getSecondBead()->_dbIndex;
        
        kstr[i] = l->getMLinker()->getStretchingConstant();
        eql[i] = l->getMLinker()->getEqLength();
        pos1[i] = l->getFirstPosition();
        pos2[i] = l->getSecondPosition();
        
        i++;
    }
}

template<class LStretchingInteractionType>
void LinkerStretching<LStretchingInteractionType>::deallocate() {
    
    delete beadSet;
    delete kstr;
    delete eql;
    delete pos1;
    delete pos2;
}


template <class LStretchingInteractionType>
double LinkerStretching<LStretchingInteractionType>::computeEnergy(double* coord, double *f, double d){
    
    double U_i;
    
    if (d == 0.0)
        U_i = _FFType.energy(coord, f, beadSet, kstr, eql, pos1, pos2);
    else
        U_i = _FFType.energy(coord, f, beadSet, kstr, eql, pos1, pos2, d);
//    std::cout<<"================="<<endl;
#ifdef CROSSCHECK
    double U2 = 0;
    double U_ii;
//    std::cout<<"NL "<<(Linker::getLinkers()).size()<<endl;
    for (auto l: Linker::getLinkers()) {
        
        Bead* b1 = l->getFirstCylinder()->getFirstBead();
        Bead* b2 = l->getFirstCylinder()->getSecondBead();
        Bead* b3 = l->getSecondCylinder()->getFirstBead();
        Bead* b4 = l->getSecondCylinder()->getSecondBead();
        double kStretch = l->getMLinker()->getStretchingConstant();
        double eqLength = l->getMLinker()->getEqLength();
        double pos1 = l->getFirstPosition();
        double pos2 = l->getSecondPosition();
        
        if (d == 0.0)
            U_ii = _FFType.energy(b1, b2, b3, b4, pos1, pos2, kStretch, eqLength);
        else
            U_ii = _FFType.energy(b1, b2, b3, b4, pos1, pos2, kStretch, eqLength, d);
        
        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_ii != U_ii || U_ii < -1.0) {
            
            U2=-1;
            break;
        }
        else
            U2 += U_ii;
    }
    if(abs(U_i-U2)<=U2/100000000000)
        std::cout<<"E L YES "<<endl;
    else
    {   std::cout<<U_i<<" "<<U2<<endl;
        exit(EXIT_FAILURE);
    }

#endif
    
    return U_i;
}

template <class LStretchingInteractionType>
void LinkerStretching<LStretchingInteractionType>::computeForces(double *coord, double *f) {
    
    _FFType.forces(coord, f, beadSet, kstr, eql, pos1, pos2);
#ifdef CROSSCHECK
    for (auto l: Linker::getLinkers()) {
        
        Bead* b1 = l->getFirstCylinder()->getFirstBead();
        Bead* b2 = l->getFirstCylinder()->getSecondBead();
        Bead* b3 = l->getSecondCylinder()->getFirstBead();
        Bead* b4 = l->getSecondCylinder()->getSecondBead();
        double kStretch = l->getMLinker()->getStretchingConstant();
        double eqLength = l->getMLinker()->getEqLength();
        
        double pos1 = l->getFirstPosition();
        double pos2 = l->getSecondPosition();
        
        if(cross_checkclass::Aux)
        {double f0 = _FFType.forcesAux(b1, b2, b3, b4, pos1, pos2, kStretch, eqLength);
//            l->getMLinker()->stretchForce = f0;
        }
        else
        {double f0 = _FFType.forces(b1, b2, b3, b4, pos1, pos2, kStretch, eqLength);
//            l->getMLinker()->stretchForce = f0;
        }
    }
    if(cross_checkclass::Aux){
        auto state=cross_check::crosscheckAuxforces(f);
        std::cout<<"F S+B+L YES "<<state<<endl;}
    else{
        auto state=cross_check::crosscheckforces(f);
        std::cout<<"F S+B+L YES "<<state<<endl;}
#endif
}


///Temlate specializations
template double LinkerStretching<LinkerStretchingHarmonic>::computeEnergy(double *coord, double *f, double d);
template void LinkerStretching<LinkerStretchingHarmonic>::computeForces(double *coord, double *f);
template void LinkerStretching<LinkerStretchingHarmonic>::vectorize();
template void LinkerStretching<LinkerStretchingHarmonic>::deallocate();

