
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

#ifndef MEDYAN_Mechanics_ForceField_Volume_TriangleCylinderVolumeFF_Hpp
#define MEDYAN_Mechanics_ForceField_Volume_TriangleCylinderVolumeFF_Hpp

#include <vector>

#include "common.h"

#include "Mechanics/ForceField/ForceField.h"

//FORWARD DECLARATIONS
class TriangleCylinderVolumeInteractions;

/// An implementation of the ForceField class that calculates
/// Triangle vs bead volume interactions.
class TriangleCylinderVolumeFF : public ForceField {
    
private:
    std::vector< std::unique_ptr< TriangleCylinderVolumeInteractions > >
        _triangleCylinderVolInteractionVector;  ///< Vector of initialized volume interactions
    
    TriangleCylinderVolumeInteractions* _culpritInteraction; ///< Culprit in case of error
public:
    /// Initialize the volume forcefields
    TriangleCylinderVolumeFF(string& interaction);

    virtual void vectorize() override { }
    virtual void cleanup() override { }

    virtual string getName() override {return "Triangle Cylinder Volume";}
    virtual void whoIsCulprit() override;

    virtual floatingpoint computeEnergy(floatingpoint* coord, bool stretched) override;
    virtual void computeForces(floatingpoint* coord, floatingpoint* f) override;

    virtual void computeLoadForces() override;
    virtual void computeLoadForce(Cylinder* c, LoadForceEnd end) const override;

    virtual vector<NeighborList*> getNeighborLists() override;
};

#endif
