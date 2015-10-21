
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

#ifndef M3SYM_BranchingFF_h
#define M3SYM_BranchingFF_h

#include <vector>

#include "common.h"

#include "ForceField.h"

//FORWARD DECLARATIONS
class BranchingInteractions;
class BranchingPoint;

/// Branching FF is an implementation of the ForceField class that
/// calculates BranchingPoint interactions.
class BranchingFF : public ForceField {
    
private:
    vector<unique_ptr<BranchingInteractions>>
    _branchingInteractionVector; ///< Vector of initialized branching interactions
    
    BranchingInteractions* _culpritInteraction; ///< Culprit in case of error
public:
    /// Constructor, intializes all interaction at the branching point
    BranchingFF(string& stretching, string& bending,
                string& dihedral, string& position);
    
    virtual string getName() {return "Branching";}
    virtual void whoIsCulprit();
    
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();
    
    virtual vector<NeighborList*> getNeighborLists() {return vector<NeighborList*>{};}
};

#endif
