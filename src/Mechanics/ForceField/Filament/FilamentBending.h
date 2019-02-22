
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
    double *kbend;
    double *eqt;

#ifdef CUDAACCL
    int * gpu_beadSet;
    double *gpu_kbend;
    double *gpu_eqt;
    int * gpu_params;
    CUDAvars cvars;
    double *F_i;
    cudaStream_t stream = NULL;
#endif
public:
    
    ///Array describing indexed set of interactions
    ///For filaments, this is a 3-bead potential
    const static int n = 3;
    
    virtual void vectorize();
    virtual void deallocate();
    
    virtual double computeEnergy(double *coord, double *f, double d);
    virtual void computeForces(double *coord, double *f);
    
    virtual const string getName() {return "Filament Bending";}
};


#endif
