
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

#include "NeighborListContainer.h"
#include "SysParams.h"
//FORWARD DECLARATIONS
class BoundaryElement;
class Bead;

/// Represents a BoundaryElement interaction with a Bead.
class BoundaryInteractions : public BBENLContainer {
    
public:
    /// Constructor, intializes the neighbor list needed
    BoundaryInteractions()
    : BBENLContainer(SysParams::Boundaries().BoundaryCutoff) {}
    
    /// Compute energy of this interaction
    virtual double computeEnergy(BoundaryElement*, Bead*, double d) = 0;
    /// Compute forces of this interaction
    virtual void computeForces(BoundaryElement*, Bead*) = 0;
    /// Compute auxiliary forces of this interaction
    virtual void computeForcesAux(BoundaryElement*, Bead*) = 0;
    
    /// Get the name of this interaction
    virtual const string getName() = 0;
};

#endif
