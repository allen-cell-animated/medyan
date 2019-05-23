
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

#ifndef MEDYAN_BubbleFF_h
#define MEDYAN_BubbleFF_h

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
    /// Initialize the forcefields
    BubbleFF(string type, string mtoc);
    
    virtual void vectorize();
    virtual void cleanup();
    
    virtual string getName() {return "Bubble";}
    virtual void whoIsCulprit();
    
    virtual double computeEnergy(double *coord, bool stretched = false) override;
    virtual void computeForces(double *coord, double *f);
//    virtual void computeForcesAux();
    
    /// BubbleFF can compute load forces on all Bead from Bubble elements
    virtual void computeLoadForces();
    
    virtual vector<NeighborList*> getNeighborLists();
};

#endif
