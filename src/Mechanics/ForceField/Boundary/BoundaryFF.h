
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

#ifndef MEDYAN_BoundaryFF_h
#define MEDYAN_BoundaryFF_h

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
    
protected:
    /// The culprit in the case of an error
    BoundaryInteractions* _culpritInteraction;
    
public:
    /// Initialize the forcefields (repulsion, attraction, etc)
    BoundaryFF(string type);
    
    virtual void vectorize();
    virtual void cleanup();

    virtual string getName() {return "Boundary";}
    virtual void whoIsCulprit();
    
    virtual double computeEnergy(double *coord, bool stretched = false) override;
    virtual void computeForces(double *coord, double *f);
    
    /// BoundaryFF can compute load forces from all boundaries.
    virtual void computeLoadForces();
    
    virtual vector<NeighborList*> getNeighborLists();
};

#endif
