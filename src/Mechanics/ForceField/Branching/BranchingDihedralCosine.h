
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

#ifndef MEDYAN_BranchingDihedralCosine_h
#define MEDYAN_BranchingDihedralCosine_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A cosine potential used by the BranchingDihedralTemplate.
class BranchingDihedralCosine {
    
public:
    double energy(double *coord, double *f, int *beadSet,
                  double *kdih, double *pos);
    
    double energy(double *coord, double *f, int *beadSet,
                  double *kdih, double *pos, double d);

    void forces(double *coord, double *f, int *beadSet,
                double *kdih, double *pos);
#ifdef CUDAACCL
    void optimalblocksnthreads(int nint);

    double* energy(double *coord, double *f, int *beadSet, double *kdih,
                   double *pos, int *params);

    double* energy(double *coord, double *f, int *beadSet, double *kdih, double *pos,
                   double *z, int *params);

    void forces(double *coord, double *f, int *beadSet, double *kdih, double *pos, int *params);
    void deallocate();
    vector<int> blocksnthreadse;
    vector<int> blocksnthreadsez;
    vector<int> blocksnthreadsf;
    vector<int> bntaddvec2;
    static void checkforculprit();
    double *gU_i;
    double *gU_sum;
    char *gFF, *ginteraction;
    cudaStream_t stream = NULL;
#endif
};

#endif
