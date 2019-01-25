
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

#ifndef MEDYAN_MTOCBendingCosine_h
#define MEDYAN_MTOCBendingCosine_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A harmonic potential used by the MTOCAttachment template.
class MTOCBendingCosine {
    
public:
    double energy(double *coord, double *f, int *beadSet,
                  double *kbend, double);
    double energy(double *coord, double *f, int *beadSet,
                  double *kbend, double, double);
    void forces(double *coord, double *f, int *beadSet,
                double *kbend, double);
    
    //void forcesAux(Bead*, Bead*, double, double);
};

#endif


