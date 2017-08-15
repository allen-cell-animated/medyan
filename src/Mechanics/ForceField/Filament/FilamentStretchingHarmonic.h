
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

#ifndef MEDYAN_FilamentStretchingHarmonic_h
#define MEDYAN_FilamentStretchingHarmonic_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A harmonic potential used by the FilamentStretching and MTOCAttachment template.
class FilamentStretchingHarmonic {
    
public:
    double energy(double *coord, double *f, int *beadSet,
                  double *kstr, double *eql);
    
    double energy(double *coord, double * f, int *beadSet,
                  double *kstr, double *eql, double d);
    
    void forces(double *coord, double *f, int *beadSet,
                double *kstr, double *eql);
};

#endif
