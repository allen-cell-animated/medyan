
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

    double energy(double *coord, double *f, int *beadSet, double *krep);
    
    double energy(double *coord, double *f, int *beadSet, double *krep, double d);
    
    void forces(double *coord, double *f, int *beadSet, double *krep);

#ifdef CUDAACCL
    double energy(double *coord, double *f, int *beadSet, double *krep, int *params, vector<int> blocksnthreads);

    double energy(double *coord, double *f, int *beadSet, double *krep, double *z, int *params, vector<int> blocksnthreads);

    void forces(double *coord, double *f, int *beadSet, double *krep, int *params, vector<int> blocksnthreads);

#endif
#ifdef CROSSCHECK
    double energy(Bead*, Bead*, Bead*, Bead*, double Krepuls);
    double energy(Bead*, Bead*, Bead*, Bead*, double Krepuls, double d);
    
    void forces(Bead*, Bead*, Bead*, Bead*, double Krepuls);
    void forcesAux(Bead*, Bead*, Bead*, Bead*, double Krepuls);
#endif
};
#endif
