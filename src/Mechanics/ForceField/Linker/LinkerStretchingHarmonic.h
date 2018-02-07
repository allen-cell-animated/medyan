
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

#ifndef MEDYAN_LinkerStretchingHarmonic_h
#define MEDYAN_LinkerStretchingHarmonic_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A harmonic potential used by the LinkerStretching template.
class LinkerStretchingHarmonic {
    
public:
    double energy(double *coord, double *f, int *beadSet,
                  double *kstr, double *eql, double *pos1, double *pos2);
    
    double energy(double *coord, double * f, int *beadSet,
                  double *kstr, double *eql, double *pos1, double *pos2, double d);
    
    void forces(double *coord, double *f, int *beadSet,
                double *kstr, double *eql, double *pos1, double *pos2);
#ifdef CUDAACCL
    void optimalblocksnthreads(int nint);

    double* energy(double *coord, double *f, int *beadSet, double *kstr, double *eql,
                   double *pos1, double *pos2, int *params);

    double* energy(double *coord, double *f, int *beadSet, double *kstr, double *eql, double *pos1, double *pos2,
                   double *z, int *params);

    void forces(double *coord, double *f, int *beadSet, double *kstr, double *eql, double *pos1, double *pos2, int
    *params);
    void deallocate();
    vector<int> blocksnthreadse;
    vector<int> blocksnthreadsez;
    vector<int> blocksnthreadsf;
    static void checkforculprit();
    double *gU_i;
    double *gU_sum;
    char *gFF, *ginteraction;
    cudaStream_t stream;
#endif
#ifdef CROSSCHECK
    double energy(Bead*, Bead*, Bead*, Bead*,
                  double position1, double position2,
                  double kStretch, double eqLength);
    double energy(Bead*, Bead*, Bead*, Bead*,
                  double position1, double position2,
                  double kStretch, double eqLength, double d);
    
    double forces(Bead*, Bead*, Bead*, Bead*,
                  double position1, double position2,
                  double kStretch, double eqLength);
    double forcesAux(Bead*, Bead*, Bead*, Bead*,
                     double position1, double position2,
                     double kStretch, double eqLength);
#endif
};

#endif
