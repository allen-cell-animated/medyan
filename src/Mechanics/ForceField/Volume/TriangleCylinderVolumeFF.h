
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

#ifndef MEDYAN_TriangleCylinderVolumeFF_h
#define MEDYAN_TriangleCylinderVolumeFF_h

#include <vector>

#include "common.h"

#include "ForceField.h"

//FORWARD DECLARATIONS
class TriangleCylinderVolumeInteractions;
class Triangle;
class Cylinder;
class Bead;

/// An implementation of the ForceField class that calculates
/// Triangle vs bead volume interactions.
class TriangleCylinderVolumeFF : public ForceField {
    
private:
    vector <unique_ptr<TriangleCylinderVolumeInteractions>>
    _triangleCylinderVolInteractionVector;  ///< Vector of initialized volume interactions
    
    TriangleCylinderVolumeInteractions* _culpritInteraction; ///< Culprit in case of error
public:
    /// Initialize the volume forcefields
    TriangleCylinderVolumeFF(string& interaction);

    virtual void vectorize() override { }
    virtual void cleanup() override { }
    
    virtual string getName() {return "Triangle Cylinder Volume";}
    virtual void whoIsCulprit();

    virtual double computeEnergy(double* coord, bool stretched) override;
    virtual void computeForces(double* coord, double* f) override;
    
    virtual void computeLoadForces();
    
    virtual vector<NeighborList*> getNeighborLists();
};

#endif
