
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef M3SYM_BoundaryCylinderRepulsion_h
#define M3SYM_BoundaryCylinderRepulsion_h

#include <vector>

#include "common.h"

#include "BoundaryInteractions.h"
#include "NeighborListImpl.h"

#include "SysParams.h"

//FORWARD DECLARATIONS
class BoundaryElement;
class Bead;

/// Represents a repulsive interaction between a BoundaryElement and Cylinder.
template <class BRepulsionInteractionType>
class BoundaryCylinderRepulsion : public BoundaryInteractions {
    
private:
    BRepulsionInteractionType _FFType;
    BoundaryCylinderNL* _neighborList; ///<Neighbor list of BoundaryElement - Cylinder
public:
    
    /// Constructor
    BoundaryCylinderRepulsion() {
        _neighborList = new BoundaryCylinderNL(SysParams::Boundaries().BoundaryCutoff);
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
    
    virtual const string getName() {return "Boundary-Cylinder Repulsion";}
};
#endif
