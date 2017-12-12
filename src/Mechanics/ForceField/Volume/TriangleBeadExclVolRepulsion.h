
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

#ifndef MEDYAN_TriangleBeadExclVolRepulsion_h
#define MEDYAN_TriangleBeadExclVolRepulsion_h

#include "common.h"

//FORWARD DECLARATIONS
class Triangle;
class Bead;

/// Represents a repulsive excluded volume potential used by the
/// CylinderExclVolume template.
class TriangleBeadExclVolRepulsion {
    
public:
    double energy(Bead*, Bead*, Bead*, Bead*, double Krepuls);
    double energy(Bead*, Bead*, Bead*, Bead*, double Krepuls, double d);
    
    void forces(Bead*, Bead*, Bead*, Bead*, double Krepuls);
    void forcesAux(Bead*, Bead*, Bead*, Bead*, double Krepuls);
};


#endif
