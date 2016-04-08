
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.0
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_BranchingStretching_h
#define MEDYAN_BranchingStretching_h

#include "common.h"

#include "BranchingInteractions.h"

//FORWARD DECLARATIONS
class BranchingPoint;

/// Represents an interaction fixing a Cylinder anchored by a BranchingPoint on the parent.
template <class BStretchingInteractionType>
class BranchingStretching : public BranchingInteractions {
    
private:
    BStretchingInteractionType _FFType;
    
public:
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();
    
    virtual const string getName() {return "Branching Stretching";}
};

#endif
