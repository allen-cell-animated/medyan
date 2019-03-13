
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

#ifndef MEDYAN_CylinderExclVolRepulsion_h
#define MEDYAN_CylinderExclVolRepulsion_h

#include "common.h"


//FORWARD DECLARATIONS
class Bead;

/// Represents a repulsive excluded volume potential used by the
/// CylinderExclVolume template.
class CylinderExclVolRepulsion {
    
public:


    floatingpoint energy(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *krep);
    
    floatingpoint energy(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *krep, floatingpoint d);
    
    void forces(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *krep);

#ifdef CUDAACCL
    void optimalblocksnthreads(int nint, cudaStream_t stream);
    void deallocate();
    floatingpoint* energy(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *krep, int *params);

    floatingpoint* energy(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *krep, floatingpoint *z, int *params);

    void forces(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *krep, int *params);
    static void checkforculprit();
    floatingpoint *gU_i;
    floatingpoint *gU_sum;
    char *gFF, *ginteraction;
    vector<int> blocksnthreadse;
    vector<int> blocksnthreadsez;
    vector<int> blocksnthreadsf;
    vector<int> bntaddvec2;
    cudaStream_t stream = NULL;

#endif
#ifdef CROSSCHECK
    floatingpoint energy(Bead*, Bead*, Bead*, Bead*, floatingpoint Krepuls);
    floatingpoint energy(Bead*, Bead*, Bead*, Bead*, floatingpoint Krepuls, floatingpoint d);
    
    void forces(Bead*, Bead*, Bead*, Bead*, floatingpoint Krepuls);
    void forcesAux(Bead*, Bead*, Bead*, Bead*, floatingpoint Krepuls);
#endif
};
#endif
