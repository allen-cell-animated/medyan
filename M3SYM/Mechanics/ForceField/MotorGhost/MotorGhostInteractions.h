
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
    
friend class MotorGhostFF;
    
protected:
    /// The motor ghost in the case of an error
    MotorGhost* _motorCulprit;
    
public:
    /// Compute the energy of this interaction
    virtual double computeEnergy(double d) = 0;
    /// Compute the forces of this interaction
    virtual void computeForces() = 0;
    /// Compute the auxiliary forces of this interaction
    virtual void computeForcesAux() = 0;

    /// Get the name of this interaction
    virtual const string getName() = 0;
};

#endif
