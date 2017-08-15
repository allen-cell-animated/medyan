
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
};

#endif
