
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

#ifndef MEDYAN_BranchingBendingCosine_h
#define MEDYAN_BranchingBendingCosine_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A cosine potential used by the BranchingBending template.
class BranchingBendingCosine {
    
public:
    double energy(double *coord, double *f, int *beadSet,
                  double *kbend, double *eqt);
    
    double energy(double *coord, double *f, int *beadSet,
                  double *kbend, double *eqt, double d);

    void forces(double *coord, double *f, int *beadSet,
                double *kbend, double *eqt);
#ifdef CUDAACCL
    void optimalblocksnthreads(int nint);

    double* energy(double *coord, double *f, int *beadSet, double *kbend, double *eqt, int *params);

    double* energy(double *coord, double *f, int *beadSet, double *kbend, double *eqt, double *z, int *params);

    void forces(double *coord, double *f, int *beadSet, double *kbend, double *eqt, int *params);
    void deallocate();
    static void checkforculprit();
    double *gU_i;
    double *gU_sum;
    char *gFF, *ginteraction;
    vector<int> blocksnthreadse;
    vector<int> blocksnthreadsez;
    vector<int> blocksnthreadsf;
    vector<int> bntaddvec2;
    cudaStream_t stream = NULL;
#endif
};

#endif
