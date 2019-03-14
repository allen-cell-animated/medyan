
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

#ifndef MEDYAN_FilamentBendingHarmonic_h
#define MEDYAN_FilamentBendingHarmonic_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A harmonic potential used by the FilamentBending template.
class FilamentBendingHarmonic {
    
public:
    double energy(double *coord, double *f, int *beadSet,
                  double *kbend, double *eqt);
    
    double energy(double *coord, double * f, int *beadSet,
                  double *kbend, double *eqt, double d);
    
    void forces(double *coord, double *f, int *beadSet,
                double *kbend, double *eqt);
#ifdef CUDAACCL
    void optimalblocksnthreads(int nint, cudaStream_t stream);
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
#ifdef CROSSCHECK
    double energy(Bead*, Bead*, Bead*, double, double);
    double energy(Bead*, Bead*, Bead*, double, double, double);
    
    void forces(Bead*, Bead*, Bead*, double, double);
    void forcesAux(Bead*, Bead*, Bead*, double, double);
#endif
};

#endif
