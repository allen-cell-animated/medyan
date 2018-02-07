
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

#ifndef MEDYAN_BranchingStretchingHarmonic_h
#define MEDYAN_BranchingStretchingHarmonic_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// Represents a harmonic potential used by the [BranchingStretching](@ref
/// BranchingStretching) template.
class BranchingStretchingHarmonic {
    
public:
    double energy(double *coord, double *f, int *beadSet,
                  double *kstr, double *eql, double *pos);
    
    double energy(double *coord, double * f, int *beadSet,
                  double *kstr, double *eql, double *pos, double d);
    
    void forces(double *coord, double *f, int *beadSet,
                double *kstr, double *eql, double *pos);
#ifdef CUDAACCL
    void optimalblocksnthreads(int nint);

    double* energy(double *coord, double *f, int *beadSet, double *kstr, double *eql, double *pos, int *params);

    double* energy(double *coord, double *f, int *beadSet, double *kstr, double *eql, double *pos, double *z, int
            *params);

    void forces(double *coord, double *f, int *beadSet, double *kstr, double *eql, double *pos, int *params);
    void deallocate();
    static void checkforculprit();
    double *gU_i;
    double *gU_sum;
    char *gFF, *ginteraction;
    vector<int> blocksnthreadse;
    vector<int> blocksnthreadsez;
    vector<int> blocksnthreadsf;
    cudaStream_t stream;
#endif
    
};

#endif
