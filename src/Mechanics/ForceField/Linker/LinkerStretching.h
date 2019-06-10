
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_LinkerStretching_h
#define MEDYAN_LinkerStretching_h

#include "common.h"
#ifdef CUDAACCL
#include "CUDAcommon.h"
#endif

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
    floatingpoint *kstr;
    floatingpoint *eql;
    floatingpoint *pos1;
    floatingpoint *pos2;
    floatingpoint *stretchforce;

#ifdef CUDAACCL
    int * gpu_beadSet;
    floatingpoint * gpu_kstr;
    floatingpoint *gpu_eql;
    int * gpu_params;
    floatingpoint *gpu_pos1;
    floatingpoint *gpu_pos2;
//    CUDAvars cvars;
    floatingpoint *F_i;
    floatingpoint *gpu_Lstretchforce;
    cudaStream_t  stream = NULL;
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
    
    
    virtual floatingpoint computeEnergy(floatingpoint *coord, floatingpoint *f, floatingpoint d);
    virtual void computeForces(floatingpoint *coord, floatingpoint *f);
    
    virtual const string getName() {return "Linker Stretching";}

    virtual void assignforcemags();
};

#endif
