
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

#ifndef MEDYAN_CylinderVolumeInteractions_h
#define MEDYAN_CylinderVolumeInteractions_h

#include "HybridNeighborListImpl.h"
#include "common.h"

//FORWARD DECLARATIONS
class NeighborList;
class HybridNeighborList;
class Cylinder;

/// Represents a volume interaction between [Cylinders](@ref Cylinder).
class CylinderVolumeInteractions {

friend class CylinderVolumeFF;

public:
    //@{
    /// The cylinder culprits in the case of an error
    static Cylinder* _cylinderCulprit1;
    static Cylinder* _cylinderCulprit2;
    //@}

    ///Vectorize the bead interactions for minimization
    virtual void vectorize() = 0;
    ///Deallocate the vectorized data
    virtual void deallocate() = 0;

    /// Compute the energy of this interaction
    virtual floatingpoint computeEnergy(floatingpoint *coord) = 0;
    /// Compute the forces of this interaction
    virtual void computeForces(floatingpoint *coord, floatingpoint *f) = 0;

    /// Get the neighbor list for this interaction
    virtual NeighborList* getNeighborList() = 0;

    /// Get the name of this interaction
    virtual const string getName() = 0;

#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
    //setter to store HNeighborList pointer.
    virtual void setHNeighborList(HybridCylinderCylinderNL* Hnl) = 0;
    //getter to obtain HNeighborList pointer.
    virtual HybridCylinderCylinderNL* getHNeighborList() = 0;
#endif

};


#endif
