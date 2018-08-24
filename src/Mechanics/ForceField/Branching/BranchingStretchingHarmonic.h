
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

#ifndef MEDYAN_BranchingStretchingHarmonic_h
#define MEDYAN_BranchingStretchingHarmonic_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// Represents a harmonic potential used by the [BranchingStretching](@ref
/// BranchingStretching) template.
class BranchingStretchingHarmonic {
    
public:
    double energy(Bead*, Bead*, Bead*,
                  double position, double kStretch, double eqLength);
    double energy(Bead*, Bead*, Bead*,
                  double position, double kStretch, double eqLength, double d);
    
    double forces(Bead*, Bead*, Bead*,
                double position, double kStretch, double eqLength);
    void forcesAux(Bead*, Bead*, Bead*,
                   double position, double kStretch, double eqLength);
    
};

#endif
