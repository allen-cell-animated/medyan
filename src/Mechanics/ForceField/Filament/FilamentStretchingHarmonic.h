
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

#ifndef MEDYAN_FilamentStretchingHarmonic_h
#define MEDYAN_FilamentStretchingHarmonic_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A harmonic potential used by the FilamentStretching and MTOCAttachment template.
class FilamentStretchingHarmonic {
    
public:
    double energy(double *coord, int *beadSet,
                  double *kstr, double *eql);
    
    double energy(double *coord, double * f, int *beadSet,
                  double *kstr, double *eql, double d);
    
    void forces(double *coord, double *f, int *beadSet,
                double *kstr, double *eql);
#ifdef CUDAACCL
    void optimalblocksnthreads(int nint, cudaStream_t stream);

    double* energy(double *coord, double *f, int *beadSet, double *kstr, double *eql, int *params);

    double* energy(double *coord, double *f, int *beadSet, double *kstr, double *eql, double *z, int *params);

    void forces(double *coord, double *f, int *beadSet, double *kstr, double *eql, int *params);
    void deallocate();
    static void checkforculprit();
    vector<int> blocksnthreadse;
    vector<int> blocksnthreadsez;
    vector<int> blocksnthreadsf;
    vector<int> bntaddvec2;
    double *gU_i;
    double *gU_sum;
    char *gFF, *ginteraction;
    cudaStream_t stream = NULL;
#endif
#ifdef CROSSCHECK
    double energy(Bead*, Bead*, double, double);
    double energy(Bead*, Bead*, double, double, double);
    
    void forces(Bead*, Bead*, double, double);
    void forcesAux(Bead*, Bead*, double, double);
#endif
};

#endif
