
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef M3SYM_MotorGhostStretchingHarmonic_h
#define M3SYM_MotorGhostStretchingHarmonic_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A harmonic potential used by the MotorGhostStretching template.
class MotorGhostStretchingHarmonic {
    
public:
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
};

#endif
