
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

#ifndef M3SYM_BoundaryFF_h
#define M3SYM_BoundaryFF_h

#include <vector>

#include "common.h"

#include "ForceField.h"

//FORWARD DECLARATIONS
class BoundaryInteractions;
class BoundaryElement;
class Bead;

/// An implementation of the ForceField class that calculates BoundaryElement
/// repulsion and attraction to [Beads](@ref Bead) in the system.
class BoundaryFF : public ForceField {
    
private:
    vector<unique_ptr<BoundaryInteractions>>
    _boundaryInteractionVector; ///< Vector of initialized boundary element interactions
    
    /// The culprit in the case of an error
    BoundaryInteractions* _culpritInteraction;
public:
    /// Initialize the forcefields (repulsion, attraction, etc)
    BoundaryFF(string type);
    
    virtual string getName() {return "Boundary";}
    virtual void whoIsCulprit();
    
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();
    
    virtual vector<NeighborList*> getNeighborLists();
};

#endif
