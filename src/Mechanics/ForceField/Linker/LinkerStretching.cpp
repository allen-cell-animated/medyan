
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

template <class LStretchingInteractionType>
double LinkerStretching<LStretchingInteractionType>::computeEnergy(double d){
    
    double U = 0;
    double U_i;
    
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
            U_i = _FFType.energy(b1, b2, b3, b4, pos1, pos2, kStretch, eqLength);
        else
            U_i = _FFType.energy(b1, b2, b3, b4, pos1, pos2, kStretch, eqLength, d);
        
        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {
            
            //set culprit and return
            _linkerCulprit = l;
            
            return -1;
        }
        else
            U += U_i;
    }
    return U;
}

template <class LStretchingInteractionType>
void LinkerStretching<LStretchingInteractionType>::computeForces() {
    
    for (auto l: Linker::getLinkers()) {
        
        Bead* b1 = l->getFirstCylinder()->getFirstBead();
        Bead* b2 = l->getFirstCylinder()->getSecondBead();
        Bead* b3 = l->getSecondCylinder()->getFirstBead();
        Bead* b4 = l->getSecondCylinder()->getSecondBead();
        double kStretch = l->getMLinker()->getStretchingConstant();
        double eqLength = l->getMLinker()->getEqLength();
        
        double pos1 = l->getFirstPosition();
        double pos2 = l->getSecondPosition();
        
        double f0 = _FFType.forces(b1, b2, b3, b4, pos1, pos2, kStretch, eqLength);
        l->getMLinker()->stretchForce = f0;
    }
    double maxF = 0;
    
    //calc max force
    for(auto b: Bead::getBeads())
        maxF = max(maxF, sqrt(b->FADotFA()));
    std::cout<<"maxF "<<getName()<<" "<<maxF<<endl;
}


template <class LStretchingInteractionType>
void LinkerStretching<LStretchingInteractionType>::computeForcesAux() {

    for (auto l: Linker::getLinkers()) {
        
        Bead* b1 = l->getFirstCylinder()->getFirstBead();
        Bead* b2 = l->getFirstCylinder()->getSecondBead();
        Bead* b3 = l->getSecondCylinder()->getFirstBead();
        Bead* b4 = l->getSecondCylinder()->getSecondBead();
        double kStretch = l->getMLinker()->getStretchingConstant();
        double eqLength = l->getMLinker()->getEqLength();
        
        double pos1 = l->getFirstPosition();
        double pos2 = l->getSecondPosition();
        
        double f0 = _FFType.forcesAux(b1, b2, b3, b4, pos1, pos2, kStretch, eqLength);
        l->getMLinker()->stretchForce = f0;
    }
    double maxF = 0;
    
    //calc max force
    for(auto b: Bead::getBeads())
        maxF = max(maxF, sqrt(b->FADotFA()));
    std::cout<<"maxF "<<getName()<<" "<<maxF<<endl;
}

///Temlate specializations
template double LinkerStretching<LinkerStretchingHarmonic>::computeEnergy(double d);
template void LinkerStretching<LinkerStretchingHarmonic>::computeForces();
template void LinkerStretching<LinkerStretchingHarmonic>::computeForcesAux();
