
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
#ifndef __M3SYM__BranchingPosition__
#define __M3SYM__BranchingPosition__

#include "common.h"

#include "BranchingInteractions.h"

//FORWARD DECLARATIONS
class BranchingPoint;

/// Represents an interaction fixing a Cylinder fixed by a BranchingPoint on the main.
template <class BStretchingInteractionType>
class BranchingPosition : public BranchingInteractions {
    
private:
    BStretchingInteractionType _FFType;
    
public:
    virtual double computeEnergy(BranchingPoint*, double d);
    virtual void computeForces(BranchingPoint*);
    virtual void computeForcesAux(BranchingPoint*);
};

#endif

