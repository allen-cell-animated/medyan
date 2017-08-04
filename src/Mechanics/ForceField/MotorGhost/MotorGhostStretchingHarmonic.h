
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

#ifndef MEDYAN_MotorGhostStretchingHarmonic_h
#define MEDYAN_MotorGhostStretchingHarmonic_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A harmonic potential used by the MotorGhostStretching template.
class MotorGhostStretchingHarmonic {
    
public:
    inline double energy(double *coord, double *f, int *beadSet,
                         double *kstr, double *eql, double *pos1, double *pos2);
    
    inline double energy(double *coord, double * f, int *beadSet,
                         double *kstr, double *eql, double *pos1, double *pos2, double d);
    
    inline void forces(double *coord, double *f, int *beadSet,
                       double *kstr, double *eql, double *pos1, double *pos2);
};

#endif
