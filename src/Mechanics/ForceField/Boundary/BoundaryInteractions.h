
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_BoundaryInteractions_h
#define MEDYAN_BoundaryInteractions_h

#include "common.h"
#include "Mechanics/ForceField/Types.hpp"

//FORWARD DECLARATIONS
class NeighborList;
class BoundaryElement;
class Bead;
class Component;
class Cylinder;

/// Represents a BoundaryElement interaction with a Bead.
class BoundaryInteractions {
    
friend class BoundaryFF;
    
public:
    using LoadForceEnd = ForceFieldTypes::LoadForceEnd;

    //@{
    /// In the case of an error
    static BoundaryElement* _boundaryElementCulprit;
    static Component* _otherCulprit;
    //@}

    virtual ~BoundaryInteractions() = default;

    ///Vectorize the bead interactions for minimization
    virtual void vectorize() = 0;
    ///Deallocate the vectorized data
    virtual void deallocate() = 0;
    
    /// Compute the energy of this interaction
    virtual floatingpoint computeEnergy(floatingpoint *coord) = 0;
    /// Compute the forces of this interaction
    virtual void computeForces(floatingpoint *coord, floatingpoint *f) = 0;
    
    /// Compute the load forces on beads from this interaction
    virtual void computeLoadForces() = 0;
    virtual void computeLoadForce(Cylinder* c, LoadForceEnd end) const { }
    
    /// Get the neighbor list for this interaction
    virtual NeighborList* getNeighborList() = 0;
    
    /// Get the name of this interaction
    virtual const string getName() = 0;
};

#endif
