
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

#ifndef MEDYAN_TriangleBeadVolumeFF_h
#define MEDYAN_TriangleBeadVolumeFF_h

#include <vector>

#include "common.h"

#include "ForceField.h"

//FORWARD DECLARATIONS
class TriangleBeadVolumeInteractions;
class Triangle;
class Bead;

/// An implementation of the ForceField class that calculates
/// Triangle vs bead volume interactions.
class TriangleBeadVolumeFF : public ForceField {
    
private:
    vector <unique_ptr<TriangleBeadVolumeInteractions>>
    _triangleBeadVolInteractionVector;  ///< Vector of initialized volume interactions
    
    TriangleBeadVolumeInteractions* _culpritInteraction; ///< Culprit in case of error
public:
    /// Initialize the volume forcefields
    TriangleBeadVolumeFF(string& interaction);
    
    virtual string getName() {return "Triangle Volume";}
    virtual void whoIsCulprit();

    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();
    
    virtual void computeLoadForces() {return;}
    
    virtual vector<NeighborList*> getNeighborLists();
};

#endif
