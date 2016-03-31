
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

#ifndef M3SYM_CylinderExclVolRepulsion_h
#define M3SYM_CylinderExclVolRepulsion_h

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// Represents a repulsive excluded volume potential used by the
/// CylinderExclVolume template.
class CylinderExclVolRepulsion {
    
public:
    double energy(Bead*, Bead*, Bead*, Bead*, double Krepuls);
    double energy(Bead*, Bead*, Bead*, Bead*, double Krepuls, double d);
    
    void forces(Bead*, Bead*, Bead*, Bead*, double Krepuls);
    void forcesAux(Bead*, Bead*, Bead*, Bead*, double Krepuls);
};


#endif
