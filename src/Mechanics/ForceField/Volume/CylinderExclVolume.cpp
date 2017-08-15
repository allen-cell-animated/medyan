
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

#include "CylinderExclVolume.h"

#include "CylinderExclVolRepulsion.h"

#include "Cylinder.h"
#include "Bead.h"

#include "MathFunctions.h"

using namespace mathfunc;

template <class CVolumeInteractionType>
int CylinderExclVolume<CVolumeInteractionType>::numInteractions;

template <class CVolumeInteractionType>
void CylinderExclVolume<CVolumeInteractionType>::vectorize() {
    
    //count interactions
    int nint = 0;
    
    for(auto ci : Cylinder::getCylinders()) {
        //do not calculate exvol for a non full length cylinder
        if(!ci->isFullLength()) continue;
        for(auto &cn : _neighborList->getNeighbors(ci))
            nint++;
    }
    
    numInteractions = nint;
    
    beadSet = new int[n * nint];
    krep = new double[nint];
    
    
    int nc = Cylinder::getCylinders().size();
    int i = 0;
    
    for (i = 0; i < nc; i++) {
        
        auto ci = Cylinder::getCylinders()[i];
        int nn = _neighborList->getNeighbors(ci).size();
        
        for (int ni = 0; ni < nn; ni++) {
            
            auto cin = _neighborList->getNeighbors(ci)[ni];
            
            beadSet[n * (i + ni)] = ci->getFirstBead()->_dbIndex;
            beadSet[n * (i + ni) + 1] = ci->getSecondBead()->_dbIndex;
            beadSet[n * (i + ni) + 2] = cin->getFirstBead()->_dbIndex;
            beadSet[n * (i + ni) + 3] = cin->getSecondBead()->_dbIndex;
            
            krep[i + ni] = ci->getMCylinder()->getExVolConst();
        }
    }
}


template <class CVolumeInteractionType>
void CylinderExclVolume<CVolumeInteractionType>::deallocate() {
    
    delete beadSet;
    delete krep;
}


template <class CVolumeInteractionType>
double CylinderExclVolume<CVolumeInteractionType>::computeEnergy(double *coord, double *f, double d) {
    
    
    double U_i = 0;
    
    if (d == 0.0)
        U_i = _FFType.energy(coord, f, beadSet, krep);
    else
        U_i = _FFType.energy(coord, f, beadSet, krep, d);
    
    return U_i;
}

template <class CVolumeInteractionType>
void CylinderExclVolume<CVolumeInteractionType>::computeForces(double *coord, double *f) {

    _FFType.forces(coord, f, beadSet, krep);
}

///Template specializations
template double CylinderExclVolume<CylinderExclVolRepulsion>::computeEnergy(double *coord, double *f, double d);
template void CylinderExclVolume<CylinderExclVolRepulsion>::computeForces(double *coord, double *f);
template void CylinderExclVolume<CylinderExclVolRepulsion>::vectorize();
template void CylinderExclVolume<CylinderExclVolRepulsion>::deallocate();


