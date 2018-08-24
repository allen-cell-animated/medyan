
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

#ifndef MEDYAN_BoundaryBubbleRepulsion_h
#define MEDYAN_BoundaryBubbleRepulsion_h

#include <vector>

#include "common.h"

#include "BoundaryInteractions.h"
#include "NeighborListImpl.h"

#include "SysParams.h"

//FORWARD DECLARATIONS
class BoundaryElement;
class Bead;

/// Represents a repulsive interaction between a BoundaryElement and Bubble.
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
   
    virtual void computeForces();
    virtual void computeForcesAux();
    
    virtual void computeLoadForces() {return;}
    
    /// Get the neighbor list for this interaction
    virtual NeighborList* getNeighborList() {return _neighborList;}
    
    virtual const string getName() {return "Boundary-Bubble Repulsion";}
};
#endif
