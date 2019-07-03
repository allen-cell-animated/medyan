
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

#ifndef MEDYAN_BoundaryCylinderRepulsion_h
#define MEDYAN_BoundaryCylinderRepulsion_h

#include <vector>

#include "common.h"
#ifdef CUDAACCL
#include "CUDAcommon.h"
#endif
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
    
    int *beadSet;
    
    ///Array describing the constants in calculation
    floatingpoint *krep;
    floatingpoint *slen;
    floatingpoint *U_i;
    int nint = 0;
    ///Array describing the number of neighbors for each boundary element (num boundary elements long)
    int *nneighbors;
#ifdef CUDAACCL
    floatingpoint *gU;
    int *gpu_beadSet;
    floatingpoint *gpu_krep;
    floatingpoint *gpu_slen;
    floatingpoint *gpu_U_i;
    int *gpu_params;
    floatingpoint *gpu_beListplane;
    int *gpu_nintperbe;
//    CUDAvars cvars;
    floatingpoint *F_i;
    cudaStream_t stream = NULL;
#endif
    
public:
    
    ///Array describing indexed set of interactions
    ///For filaments, this is a 1-bead potential
    const static int n = 1;
    
    /// Constructor
    BoundaryCylinderRepulsion() {
        _neighborList = new BoundaryCylinderNL(SysParams::Boundaries().BoundaryCutoff);
    }
    
    virtual void vectorize();
    virtual void deallocate();
    
    virtual floatingpoint computeEnergy(floatingpoint *coord, floatingpoint *f, floatingpoint d);
    //@{
    /// This repulsive force calculation also updates load forces
    /// on beads within the interaction range.
    virtual void computeForces(floatingpoint *coord, floatingpoint *f);
    
    virtual void computeLoadForces();
    //@}
    
    /// Get the neighbor list for this interaction
    virtual NeighborList* getNeighborList() {return _neighborList;}
    
    virtual const string getName() {return "Boundary-Cylinder Repulsion";}
};
#endif
