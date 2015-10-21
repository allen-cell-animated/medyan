
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

#ifndef M3SYM_BoundaryInteractions_h
#define M3SYM_BoundaryInteractions_h

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
    BoundaryElement* _boundaryElementCulprit;
    Component* _otherCulprit;
    //@}
    
public:
    /// Compute energy of this interaction
    virtual double computeEnergy(double d) = 0;
    /// Compute forces of this interaction
    virtual void computeForces() = 0;
    /// Compute auxiliary forces of this interaction
    virtual void computeForcesAux() = 0;
    
    /// Get the neighbor list for this interaction
    virtual NeighborList* getNeighborList() = 0;
    
    /// Get the name of this interaction
    virtual const string getName() = 0;
};

#endif
