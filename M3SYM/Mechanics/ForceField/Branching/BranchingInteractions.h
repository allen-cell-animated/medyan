
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

#ifndef M3SYM_BranchingInteractions_h
#define M3SYM_BranchingInteractions_h

#include <iostream>

#include "common.h"

//FORWARD DECLARATIONS
class BranchingPoint;

/// Represents an interaction at a BranchingPoint
class BranchingInteractions {
    
public:
    /// Compute the energy of this interaction
    virtual double computeEnergy(BranchingPoint*, double d) = 0;
    /// Compute forces of this interaction
    virtual void computeForces(BranchingPoint*) = 0;
    /// Compute auxiliary forces of this interaction
    virtual void computeForcesAux(BranchingPoint*) = 0;
    
    /// Get the name of this interaction
    virtual const string getName() = 0;
};

#endif
