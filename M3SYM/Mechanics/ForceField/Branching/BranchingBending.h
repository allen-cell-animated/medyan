
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

#ifndef M3SYM_BranchingBending_h
#define M3SYM_BranchingBending_h

#include "common.h"

#include "BranchingInteractions.h"

//FORWARD DECLARATIONS
class BranchingPoint;

/// Represents an interaction maintaining a BranchingPoint angle (~270 for Arp2/3)
template <class BBendingInteractionType>
class BranchingBending : public BranchingInteractions {
    
private:
    BBendingInteractionType _FFType;
    
public:
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();
    
    virtual const string getName() {return "Branching Bending";}
};

#endif
