
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2
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

//FORWARD DECLARATIONS
class NeighborList;
class BoundaryElement;
class Bead;
class Component;

/// Represents a BoundaryElement interaction with a Bead.
class BoundaryInteractions {
    
friend class BoundaryFF;
    
protected:
    //@{
    /// In the case of an error
    BoundaryElement* _boundaryElementCulprit = nullptr;
    Component* _otherCulprit = nullptr;
    //@}
    
public:
    /// Compute energy of this interaction
    virtual double computeEnergy(double d) = 0;
    /// Compute forces of this interaction
    virtual void computeForces() = 0;
    /// Compute auxiliary forces of this interaction
    virtual void computeForcesAux() = 0;
    
    /// Compute the load forces on beads from this interaction
    virtual void computeLoadForces() = 0;
    
    /// Get the neighbor list for this interaction
    virtual NeighborList* getNeighborList() = 0;
    
    /// Get the name of this interaction
    virtual const string getName() = 0;
};

#endif
