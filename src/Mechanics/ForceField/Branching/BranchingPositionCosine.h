
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

#ifndef MEDYAN_BranchingPositionCosine_h
#define MEDYAN_BranchingPositionCosine_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A cosine potential used by the BranchingPosition template.
class BranchingPositionCosine {
    
public:
    floatingpoint energy(floatingpoint *coord, int *beadSet,
                  floatingpoint *kpos, floatingpoint *pos);
    
    [[deprecated]] floatingpoint energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                  floatingpoint *kpos, floatingpoint *pos, floatingpoint d);
    
    void forces(floatingpoint *coord, floatingpoint *f, int *beadSet,
                floatingpoint *kpos, floatingpoint *pos);
#ifdef CUDAACCL
    void optimalblocksnthreads(int nint);
    floatingpoint* energy(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *kpos, floatingpoint *pos, int *params);

    floatingpoint* energy(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *kpos, floatingpoint *pos, floatingpoint *z, int *params);

    void forces(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *kpos, floatingpoint *pos, int *params);
    void deallocate();
    vector<int> blocksnthreadse;
    vector<int> blocksnthreadsez;
    vector<int> blocksnthreadsf;
    vector<int> bntaddvec2;
    static void checkforculprit();
    floatingpoint *gU_i;
    floatingpoint *gU_sum;
    char *gFF, *ginteraction;
    cudaStream_t stream = NULL;
#endif
};

#endif
