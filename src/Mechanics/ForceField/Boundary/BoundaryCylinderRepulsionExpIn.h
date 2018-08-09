
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

#ifndef MEDYAN_BoundaryCylinderRepulsionExpIn_h
#define MEDYAN_BoundaryCylinderRepulsionExpIn_h

#include <vector>
#include <cmath>

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A exponential repulsive potential used by the BoundaryCylinderRepulsion template.
class BoundaryCylinderRepulsionExpIn {
    
    public:
    double energy(double *coord, double *f, int *beadSet,
                  double *krep, double *slen, int *nneighbors);
    
    double energy(double *coord, double *f, int *beadSet,
                  double *krep, double *slen, int *nnneighbors, double d);
    
    void forces(double *coord, double *f, int *beadSet,
                double *krep, double *slen, int *nneighbors);
    
    double loadForces(double r, double krep , double slen);
    
#ifdef CUDAACCL
    void optimalblocksnthreads(int nint);
    
    double* energy(double *coord, double *f, int *beadSet, double *krep, double *slen,
                   int* nintvec, double* beListplane, int *params);
    
    double* energy(double *coord, double *f, int *beadSet, double *krep, double *slen,
                   int* nintvec, double* beListplane, double *z, int *params);
    
    void forces(double *coord, double *f, int *beadSet, double *krep, double *slen,
                int* nintvec, double* beListplane, int *params);
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
    private:
#ifdef CUDAACCL
    //    double *F_i;
    //    double *forcecopy;
#endif
};

#endif

