
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

#ifndef MEDYAN_LinkerStretching_h
#define MEDYAN_LinkerStretching_h

#include "common.h"
#include "CUDAcommon.h"

#include "LinkerInteractions.h"
#include "Linker.h"

//FORWARD DECLARATIONS
class Linker;

/// Represents a Linker stretching interaction
template <class LStretchingInteractionType>
class LinkerStretching : public LinkerInteractions {
    
private:
    LStretchingInteractionType _FFType;

    int *beadSet;
    
    ///Array describing the constants in calculation
    double *kstr;
    double *eql;
    double *pos1;
    double *pos2;

#ifdef CUDAACCL
    int * gpu_beadSet;
    double * gpu_kstr;
    double *gpu_eql;
    int * gpu_params;
    double *gpu_pos1;
    double *gpu_pos2;
//    CUDAvars cvars;
    double *F_i;
#endif
    
public:
    
    ///Array describing indexed set of interactions
    ///For linkers, this is a 4-bead potential
    const static int n = 4;
    
    ///< Constructor
    LinkerStretching () {}
    ~LinkerStretching () {}

    virtual void vectorize();
    virtual void deallocate();
    
    
    virtual double computeEnergy(double *coord, double *f, double d);
    virtual void computeForces(double *coord, double *f);
    
    virtual const string getName() {return "Linker Stretching";}
};

#endif
