
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_FilamentBending_h
#define MEDYAN_FilamentBending_h

#include "common.h"
#ifdef CUDAACCL
#include "CUDAcommon.h"
#endif
#include "FilamentInteractions.h"
//FORWARD DECLARATIONS
class Filament;

/// Represents a Filament bending interaction
template <class FBendingInteractionType>
class FilamentBending : public FilamentInteractions {
    
private:
    FBendingInteractionType _FFType;
    
    int *beadSet;
    
    ///Array describing the constants in calculation
    floatingpoint *kbend;
    floatingpoint *eqt;

#ifdef CUDAACCL
    int * gpu_beadSet;
    floatingpoint *gpu_kbend;
    floatingpoint *gpu_eqt;
    int * gpu_params;
    CUDAvars cvars;
    floatingpoint *F_i;
    cudaStream_t stream = NULL;
#endif
public:
    
    ///Array describing indexed set of interactions
    ///For filaments, this is a 3-bead potential
    const static int n = 3;
    
    virtual void vectorize();
    virtual void deallocate();
    
    virtual floatingpoint computeEnergy(floatingpoint *coord, floatingpoint *f, floatingpoint d);
    virtual void computeForces(floatingpoint *coord, floatingpoint *f);
    
    virtual const string getName() {return "Filament Bending";}
};


#endif
