
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

#ifndef M3SYM_BoundaryRepulsion_h
#define M3SYM_BoundaryRepulsion_h

#include <vector>

#include "common.h"

#include "BoundaryInteractions.h"

//FORWARD DECLARATIONS
class BoundaryElement;
class Bead;

/// Represents a repulsive interaction between a BoundaryElement and Bead.
template <class BRepulsionInteractionType>
class BoundaryRepulsion : public BoundaryInteractions {
    
private:
    BRepulsionInteractionType _FFType;
    
public:
    virtual double computeEnergy(BoundaryElement*, Bead*, double d);
    //@{
    /// This repulsive force calculation also updates load forces
    /// on beads within the interaction range.
    virtual void computeForces(BoundaryElement*, Bead*);
    virtual void computeForcesAux(BoundaryElement*, Bead*);
    //@}
    
    virtual const string getName() {return "Repulsion";}
};
#endif
