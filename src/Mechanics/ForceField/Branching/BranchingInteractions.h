
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_BranchingInteractions_h
#define MEDYAN_BranchingInteractions_h

#include <iostream>

#include "common.h"

//FORWARD DECLARATIONS
class BranchingPoint;

/// Represents an interaction at a BranchingPoint
class BranchingInteractions {
    
friend class BranchingFF;
    
protected:
    /// The branching point in the case of an error
    BranchingPoint* _branchingCulprit = nullptr;
    
public:
    /// Compute the energy of this interaction
    virtual double computeEnergy(double d) = 0;
    /// Compute forces of this interaction
    virtual void computeForces() = 0;
    /// Compute auxiliary forces of this interaction
    virtual void computeForcesAux() = 0;
    
    /// Get the name of this interaction
    virtual const string getName() = 0;
};

#endif
