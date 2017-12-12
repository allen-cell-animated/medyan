
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2017-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_TriangleCylinderBeadExclVolRepulsion_h
#define MEDYAN_TriangleCylinderBeadExclVolRepulsion_h

#include "common.h"

//FORWARD DECLARATIONS
class Triangle;
class Bead;

/// Represents a repulsive excluded volume potential used by the
/// TriangleCylinderExclVolume template.
class TriangleCylinderBeadExclVolRepulsion {
    
public:
    double energy(Triangle*, Bead*, double kExVol);
    double energy(Triangle*, Bead*, double kExVol, double d);
    
    void forces(Triangle*, Bead*, double kExVol);
    void forcesAux(Triangle*, Bead*, double kExVol);
};


#endif
