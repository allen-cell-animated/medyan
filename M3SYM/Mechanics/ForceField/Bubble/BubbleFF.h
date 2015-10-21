
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

#ifndef M3SYM_BubbleFF_h
#define M3SYM_BubbleFF_h

#include <vector>

#include "common.h"

#include "ForceField.h"

//FORWARD DECLARATIONS
class BubbleInteractions;
class Bead;

/// An implementation of the ForceField class that calculates Bubble
/// repulsion and attraction to [Beads](@ref Bead) in the system.
class BubbleFF : public ForceField {
    
private:
    vector<unique_ptr<BubbleInteractions>>
    _bubbleInteractionVector; ///< Vector of initialized bubble interactions
    
    /// The culprit in the case of an error
    BubbleInteractions* _culpritInteraction;
public:
    /// Initialize the forcefields (repulsion, attraction, etc)
    BubbleFF(string type);
    
    virtual string getName() {return "Bubble";}
    virtual void whoIsCulprit();
    
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();
    
    virtual vector<NeighborList*> getNeighborLists();
};

#endif
