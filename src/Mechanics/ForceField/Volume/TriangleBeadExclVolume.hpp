
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

#ifndef MEDYAN_Mechanics_ForceField_Volume_TriangleBeadExclVolume_Hpp
#define MEDYAN_Mechanics_ForceField_Volume_TriangleBeadExclVolume_Hpp

#include <memory> // unique_ptr

#include "Mechanics/ForceField/Volume/TriangleBeadVolumeInteractions.hpp"
#include "Structure/NeighborListImpl.h"
#include "SysParams.h"

//FORWARD DECLARATIONS
class Triangle;
class Cylinder;
class Bead;

/// Represents an excuded volume interaction between a triangle and a cylinder (bead).
template < typename InteractionType >
class TriangleBeadExclVolume : public TriangleBeadVolumeInteractions {
    
private:
    InteractionType _FFType;
    std::unique_ptr<TriangleFilBeadNL> _neighborList;  ///< Neighbor list of triangle-bead

public:
    ///Constructor
    TriangleBeadExclVolume() :
        _neighborList(std::make_unique< TriangleFilBeadNL >(
            SysParams::Mechanics().MemBeadVolumeCutoff,
            SysParams::Mechanics().MemBeadVolumeCutoffMech
        ))
    { }

    virtual void vectorize() override {}

    virtual floatingpoint computeEnergy(const floatingpoint* coord, bool stretched) override;
    virtual void computeForces(const floatingpoint* coord, floatingpoint* force) override;

    virtual void computeLoadForces() const override;
    virtual void computeLoadForce(Cylinder* c, LoadForceEnd end) const override;

    /// Get the neighbor list for this interaction
    virtual NeighborList* getNeighborList() override { return _neighborList.get(); }
    
    virtual const string getName() override {return "Triangle Bead Excluded Volume";}
};

#endif
