
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

#ifndef MEDYAN_FilamentBendingCosine_h
#define MEDYAN_FilamentBendingCosine_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A cosine potential used by the [FilamentBending](@ref FilamentBending) template.
class FilamentBendingCosine {
    
public:
    inline double energy(double *coord, double *f, int *beadSet,
                         double *kstr, double *eqt);
    
    inline double energy(double *coord, double * f, int *beadSet,
                         double *kstr, double *eqt, double d);
    
    inline void forces(double *coord, double *f, int *beadSet,
                       double *kstr, double *eqt);
};

#endif
