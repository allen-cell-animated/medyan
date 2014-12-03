
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_MotorGhostStretchingHarmonic_h
#define M3SYM_MotorGhostStretchingHarmonic_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A harmonic potential used by the [MotorStretching](@ref MotorGhostStretching) template.
class MotorGhostStretchingHarmonic {
    
public:
    double energy(Bead*, Bead*, Bead*, Bead*, double position1, double position2, double kStr, double L);
    double energy(Bead*, Bead*, Bead*, Bead*, double position1, double position2, double kStr, double L, double d);
    void forces(Bead*, Bead*, Bead*, Bead*, double position1, double position2, double kStr, double L);
    void forcesAux(Bead*, Bead*, Bead*, Bead*, double position1, double position2, double kStr, double L);
    
};

#endif
