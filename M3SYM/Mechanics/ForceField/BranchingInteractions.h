//
//  BranchingInteraction.h
//  M3SYM
//
//  Created by Konstantin Popov on 12/2/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef M3SYM_BranchingInteractions_h
#define M3SYM_BranchingInteractions_h


#include <iostream>

#include "common.h"

//FORWARD DECLARATIONS
class BranchingPoint;

/// BranchingInteractions class represents an interaction in a branching point
class BranchingInteractions {
private:
    string _name;
    
public:
    /// Compute the energy of this interaction
    virtual double computeEnergy(BranchingPoint*,  double d) = 0;
    /// Compute forces of this interaction
    virtual void computeForces(BranchingPoint*) = 0;
    /// Compute auxiliary forces of this interaction
    virtual void computeForcesAux(BranchingPoint*) = 0;
    
    /// Get the name of this interaction
    const string& getName() {return _name;}
    
};



#endif
