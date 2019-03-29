
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

#ifndef MEDYAN_MTOCAttachmentHarmonic_h
#define MEDYAN_MTOCAttachmentHarmonic_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A harmonic potential used by the MTOCAttachment template.
class MTOCAttachmentHarmonic {
    
public:
    double energy(double *coord, double *f, int *beadSet,
                  double *kstr, double);
    double energy(double *coord, double *f, int *beadSet,
                  double *kstr, double, double);
    void forces(double *coord, double *f, int *beadSet,
                double *kstr, double);
    
    //void forcesAux(Bead*, Bead*, double, double);
};

#endif

