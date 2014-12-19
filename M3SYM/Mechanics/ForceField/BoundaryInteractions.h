
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_BoundaryInteractions_h
#define M3SYM_BoundaryInteractions_h

#include "common.h"

#include "NeighborListContainer.h"
#include "SystemParameters.h"

//FORWARD DECLARATIONS
class BoundaryElement;
class Bead;

/// Represents a BoundaryElement interaction with a Bead.
class BoundaryInteractions : public BoundaryElementNLContainer {
private:
    string _name; ///< Name of interaction
    
public:
    /// Constructor, intializes the neighbor list needed
    BoundaryInteractions()
        : BoundaryElementNLContainer(SystemParameters::Boundaries().boundaryCutoff) {}
    
    /// Compute energy of this interaction
    virtual double computeEnergy(BoundaryElement*, Bead*, double d) = 0;
    /// Compute forces of this interaction
    virtual void computeForces(BoundaryElement*, Bead*) = 0;
    /// Compute auxiliary forces of this interaction
    virtual void computeForcesAux(BoundaryElement*, Bead*) = 0;
    
    /// Get name of interaction
    string getName() {return _name;}
    
};

#endif
