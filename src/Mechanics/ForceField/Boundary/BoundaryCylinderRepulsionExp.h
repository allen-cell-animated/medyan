
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

#ifndef MEDYAN_BoundaryCylinderRepulsionExp_h
#define MEDYAN_BoundaryCylinderRepulsionExp_h

#include <vector>
#include <cmath>

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A exponential repulsive potential used by the BoundaryCylinderRepulsion template.
class BoundaryCylinderRepulsionExp {
    
public:
    floatingpoint energy(floatingpoint *coord, totalforcefloatingpoint *f, int *beadSet,
                  floatingpoint *krep, floatingpoint *slen, int *nneighbors);
    
    floatingpoint energy(floatingpoint *coord, totalforcefloatingpoint *f, int *beadSet,
                  floatingpoint *krep, floatingpoint *slen, int *nnneighbors, floatingpoint d);
    
    void forces(floatingpoint *coord, totalforcefloatingpoint *f, int *beadSet,
                floatingpoint *krep, floatingpoint *slen, int *nneighbors);
    
    floatingpoint loadForces(floatingpoint r, floatingpoint krep , floatingpoint slen);

#ifdef CUDAACCL
    void optimalblocksnthreads(int nint, cudaStream_t stream);

    floatingpoint* energy(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *krep, floatingpoint *slen,
                   int* nintvec, floatingpoint* beListplane, int *params);

    floatingpoint* energy(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *krep, floatingpoint *slen,
                   int* nintvec, floatingpoint* beListplane, floatingpoint *z, int *params);

    void forces(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *krep, floatingpoint *slen,
                int* nintvec, floatingpoint* beListplane, int *params);
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
private:
#ifdef CUDAACCL
//    floatingpoint *F_i;
//    floatingpoint *forcecopy;
#endif
};

#endif
