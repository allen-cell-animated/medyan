
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


#ifndef MEDYAN_MotorGhostInteractions_h
#define MEDYAN_MotorGhostInteractions_h

#include "common.h"

//FORWARD DECLARATIONS
class MotorGhost;

/// Represents an internal MotorGhost interaction
class MotorGhostInteractions {
    
friend class MotorGhostFF;
    
protected:
    /// The motor ghost in the case of an error
    MotorGhost* _motorCulprit = nullptr;
    
public:
    /// Compute the energy of this interaction
    virtual double computeEnergy(bool stretched) = 0;
    /// Compute the forces of this interaction
    virtual void computeForces() = 0;
    /// Compute the auxiliary forces of this interaction
    virtual void computeForcesAux() = 0;

    /// Get the name of this interaction
    virtual const string getName() = 0;
};

#endif
