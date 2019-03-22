
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

#ifndef MEDYAN_BranchingDihedral_h
#define MEDYAN_BranchingDihedral_h

#include "common.h"
#include "CUDAcommon.h"
#include "BranchingInteractions.h"

//FORWARD DECLARATIONS
class BranchingPoint;

/// Represents an interaction keeping BranchingPoint in dihedral plane
template <class BDihedralInteractionType>
class BranchingDihedral : public BranchingInteractions {
    
private:
    BDihedralInteractionType _FFType;
    
    int *beadSet;
    
    ///Array describing the constants in calculation
    floatingpoint *kdih;
    floatingpoint *pos;
#ifdef CUDAACCL
    int * gpu_beadSet;
    floatingpoint * gpu_kdih;
    floatingpoint *gpu_pos;
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
    
    virtual floatingpoint computeEnergy(floatingpoint *coord, totalforcefloatingpoint *f, floatingpoint d);
    virtual void computeForces(floatingpoint *coord, totalforcefloatingpoint *f);
    
    virtual const string getName() {return "Branching Dihedral";}
};

#endif
