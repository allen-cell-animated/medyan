
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

    int *beadSet;

    ///Array describing the constants in calculation
    floatingpoint *krep;
	floatingpoint *slen;
	floatingpoint *U_i;
    int nint = 0;
    ///Array describing the number of neighbors for each boundary element (num boundary elements long)
    int *nneighbors;
public:
    
    const static int n = 1;

    /// Constructor
    BoundaryBubbleRepulsion() {
        _neighborList = new BoundaryBubbleNL(SysParams::Boundaries().BoundaryCutoff);
    }
    
    virtual floatingpoint computeEnergy(floatingpoint *coord, floatingpoint *f, floatingpoint d);
   
    virtual void computeForces(floatingpoint *coord, floatingpoint *f);
    //virtual void computeForcesAux();
    
    virtual void computeLoadForces()  {return;}
    
    /// Get the neighbor list for this interaction
    virtual NeighborList* getNeighborList() {return _neighborList;}
    
    virtual const string getName() {return "Boundary-Bubble Repulsion";}

    //TODO needs implmenetation @{
    virtual void vectorize();
    virtual void deallocate();

    //virtual double computeEnergy(double *coord, double *f, double d)override { return 0.0; }
    //@{
    /// This repulsive force calculation also updates load forces
    /// on beads within the interaction range.
    //virtual void computeForces(double *coord, double *f){};
    //@}
};
#endif
