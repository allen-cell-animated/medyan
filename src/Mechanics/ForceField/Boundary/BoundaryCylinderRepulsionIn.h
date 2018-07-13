
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_BoundaryCylinderRepulsionIn_h
#define MEDYAN_BoundaryCylinderRepulsionIn_h

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
class BoundaryCylinderRepulsionIn : public BoundaryInteractions {
    
    private:
    BRepulsionInteractionType _FFType;
    BoundaryCylinderNL* _neighborList; ///<Neighbor list of BoundaryElement - Cylinder
    
    int *beadSet;
    
    ///Array describing the constants in calculation
    double *krep;
    double *slen;
    double *U_i;
    int nint = 0;
    ///Array describing the number of neighbors for each boundary element (num boundary elements long)
    int *nneighbors;
#ifdef CUDAACCL
    double *gU;
    cudaStream_t  stream;
    int *gpu_beadSet;
    double *gpu_krep;
    double *gpu_slen;
    double *gpu_U_i;
    int *gpu_params;
    double *gpu_beListplane;
    int *gpu_nintperbe;
    //    CUDAvars cvars;
    double *F_i;
#endif
    
    public:
    
    ///Array describing indexed set of interactions
    ///For filaments, this is a 1-bead potential
    const static int n = 1;
    
    /// Constructor
    BoundaryCylinderRepulsionIn() {
        _neighborList = new BoundaryCylinderNL(SysParams::Boundaries().BoundaryCutoff);
    }
    
    virtual void vectorize();
    virtual void deallocate();
    
    virtual double computeEnergy(double *coord, double *f, double d);
    //@{
    /// This repulsive force calculation also updates load forces
    /// on beads within the interaction range.
    virtual void computeForces(double *coord, double *f);
    
    virtual void computeLoadForces();
    //@}
    
    /// Get the neighbor list for this interaction
    virtual NeighborList* getNeighborList() {return _neighborList;}
    
    virtual const string getName() {return "Boundary-Cylinder Repulsion";}
};
#endif

