
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

#ifndef MEDYAN_VolumeConservationFF_h
#define MEDYAN_VolumeConservationFF_h

#include <vector>

#include "common.h"

#include "ForceField.h"

//FORWARD DECLARATIONS
class VolumeConservationInteractions;
//class Cylinder;

/// An implementation of the ForceField class that calculates
/// Cylinder volume interactions.
class VolumeConservationFF : public ForceField {
    
private:
    vector<unique_ptr<VolumeConservationInteractions>>
        _volumeConservationInteractionVector;  ///< Vector of initialized volume interactions
    
    VolumeConservationInteractions* _culpritInteraction; ///< Culprit in case of error
public:
    /// Initialize the volume forcefields
    VolumeConservationFF(string& interaction);

    virtual void vectorize() override {}
    virtual void cleanup() override {}

    virtual string getName() { return "Volume conservation"; }
    virtual void whoIsCulprit();

    virtual double computeEnergy(double* coord, bool stretched) override;
    virtual void computeForces(double* coord, double* f) override;
    
    virtual void computeLoadForces() {}
    
    virtual vector<NeighborList*> getNeighborLists() { return vector<NeighborList*>{}; }
};

#endif
