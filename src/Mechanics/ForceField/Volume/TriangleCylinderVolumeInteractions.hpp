
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2017-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_Mechanics_ForceField_Volume_TriangleCylinderVolumeInteractions_Hpp
#define MEDYAN_Mechanics_ForceField_Volume_TriangleCylinderVolumeInteractions_Hpp

#include "common.h" // floatingpoint
#include "Mechanics/ForceField/Types.hpp"

//FORWARD DECLARATIONS
class NeighborList;
class Triangle;
class Cylinder;
class Bead;

/// Represents a volume interaction between [Cylinders](@ref Cylinder).
class TriangleCylinderVolumeInteractions {

friend class TriangleCylinderVolumeFF;
    
protected:
    //@{
    /// The triangle and bead culprits in the case of an error
    Triangle* _triangleCulprit = nullptr;
    Cylinder* _cylinderCulprit = nullptr;
    //@}
    
public:
    using LoadForceEnd = ForceFieldTypes::LoadForceEnd;

    virtual ~TriangleCylinderVolumeInteractions() = default;

    /// Compute the energy of this interaction
    virtual floatingpoint computeEnergy(const floatingpoint* coord, bool stretched) = 0;
    /// Compute the forces of this interaction
    virtual void computeForces(const floatingpoint* coord, floatingpoint* force) = 0;

    /// Compute the load forces on beads from this interaction
    virtual void computeLoadForces() const = 0;
    virtual void computeLoadForce(Cylinder* c, LoadForceEnd end) const { }
    
    /// Get the neighbor list for this interaction
    virtual NeighborList* getNeighborList() = 0;
    
    /// Get the name of this interaction
    virtual const string getName() = 0;
};


#endif
