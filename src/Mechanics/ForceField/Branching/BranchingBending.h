
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

#ifndef MEDYAN_BranchingBending_h
#define MEDYAN_BranchingBending_h

#include "common.h"
#include "CUDAcommon.h"
#include "BranchingInteractions.h"

//FORWARD DECLARATIONS
class BranchingPoint;

/// Represents an interaction maintaining a BranchingPoint angle (~70 for Arp2/3)
template <class BBendingInteractionType>
class BranchingBending : public BranchingInteractions {
    
private:
    BBendingInteractionType _FFType;
    
    int *beadSet;
    
    ///Array describing the constants in calculation
    floatingpoint *kbend;
    floatingpoint *eqt;
#ifdef CUDAACCL
    int * gpu_beadSet;
    floatingpoint * gpu_kbend;
    floatingpoint *gpu_eqt;
    int * gpu_params;
    CUDAvars cvars;
    floatingpoint *F_i;
#endif
public:
    
    ///Array describing indexed set of interactions
    ///this is a 4-bead potential
    const static int n = 4;
    
    virtual void vectorize();
    virtual void deallocate();
    
    virtual floatingpoint computeEnergy(floatingpoint *coord, floatingpoint *f, floatingpoint d);
    virtual void computeForces(floatingpoint *coord, floatingpoint *f);
    
    virtual const string getName() {return "Branching Bending";}
};

#endif
