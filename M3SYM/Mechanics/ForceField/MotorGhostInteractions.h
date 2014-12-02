
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


#ifndef M3SYM_MotorGhostInteractions_h
#define M3SYM_MotorGhostInteractions_h

#include "common.h"

//FORWARD DECLARATIONS
class MotorGhost;

/// Represents an internal motor interaction
class MotorGhostInteractions {
private:
    string _name; ///< Name of interaction
    
public:
    /// Compute the energy of this interaction
    virtual double computeEnergy( MotorGhost*,  double d) = 0;
    /// Compute the forces of this interaction
    virtual void computeForces(MotorGhost*) = 0;
    /// Compute the auxiliary forces of this interaction
    virtual void computeForcesAux(MotorGhost*) = 0;

    /// Get name of this interaction
    const string& getName() {return _name;}
};

#endif
