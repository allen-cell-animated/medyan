
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------


#ifndef M3SYM_MotorGhostInteractions_h
#define M3SYM_MotorGhostInteractions_h

#include "common.h"

//FORWARD DECLARATIONS
class MotorGhost;

/// Represents an internal MotorGhost interaction
class MotorGhostInteractions {
    
public:
    /// Compute the energy of this interaction
    virtual double computeEnergy( MotorGhost*,  double d) = 0;
    /// Compute the forces of this interaction
    virtual void computeForces(MotorGhost*) = 0;
    /// Compute the auxiliary forces of this interaction
    virtual void computeForcesAux(MotorGhost*) = 0;

    /// Get the name of this interaction
    virtual const string getName() = 0;
};

#endif
