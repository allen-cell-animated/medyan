
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

#ifndef MEDYAN_FilamentStretching_h
#define MEDYAN_FilamentStretching_h
#include "Filament.h"
#include "Cylinder.h"
#include "common.h"
#include "CUDAcommon.h"
#include "FilamentInteractions.h"
/// Represents a Filament stretching interaction
template <class FStretchingInteractionType>
class FilamentStretching : public FilamentInteractions {
    
private:
    FStretchingInteractionType _FFType; 
    
    int *beadSet;
    ///Array describing the constants in calculation
    double *kstr;
    double *eql;

#ifdef CUDAACCL
    int * gpu_beadSet;
    double * gpu_kstr;
    double *gpu_eql;
    int * gpu_params;
    CUDAvars cvars;
    double *F_i;
#endif


public:
    
    ///Array describing indexed set of interactions
    ///For filaments, this is a 2-bead potential
    const static int n = 2;
    
    virtual void vectorize();
    virtual void deallocate();
    
    virtual double computeEnergy(double *coord, double *f, double d);
    virtual void computeForces(double *coord, double *f);
    
    virtual const string getName() {return "Filament Stretching";}

//    virtual void whoisCulprit();
};
#endif
