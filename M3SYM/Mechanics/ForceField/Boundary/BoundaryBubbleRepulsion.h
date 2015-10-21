
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

#ifndef M3SYM_BoundaryBubbleRepulsion_h
#define M3SYM_BoundaryBubbleRepulsion_h

#include <vector>

#include "common.h"

#include "BoundaryInteractions.h"
#include "NeighborListImpl.h"

#include "SysParams.h"

//FORWARD DECLARATIONS
class BoundaryElement;
class Bead;

/// Represents a repulsive interaction between a BoundaryElement and Bubbles.
template <class BRepulsionInteractionType>
class BoundaryBubbleRepulsion : public BoundaryInteractions {
    
private:
    BRepulsionInteractionType _FFType;
    BoundaryBubbleNL* _neighborList; ///<Neighbor list of Bubble's bead - BoundaryElement
public:
    
    /// Constructor
    BoundaryBubbleRepulsion() {
        _neighborList = new BoundaryBubbleNL(SysParams::Boundaries().BoundaryCutoff);
    }
    
    virtual double computeEnergy(double d);
    //@{
    /// This repulsive force calculation also updates load forces
    /// on beads within the interaction range.
    virtual void computeForces();
    virtual void computeForcesAux();
    //@}
    
    /// Get the neighbor list for this interaction
    virtual NeighborList* getNeighborList() {return _neighborList;}
    
    virtual const string getName() {return "Boundary-Bubble Repulsion";}
};
#endif
