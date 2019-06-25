
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

#ifndef MEDYAN_MembraneFF_h
#define MEDYAN_MembraneFF_h

#include "common.h"

#include "Mechanics/ForceField/ForceField.h"

//FORWARD DECLARATIONS
class MembraneInteractions;
class Membrane;

/// An implementation of the ForceField class that calculates membrane interactions.
class MembraneFF : public ForceField {
 
private:
    std::vector< std::unique_ptr< MembraneInteractions > >
        _membraneInteractionVector; ///< Vector of initialized membrane interactions
    
    MembraneInteractions* _culpritInteraction; ///< Culprit in case of error
public:
    /// Constructor, intializes stretching and bending forces
    MembraneFF(const string& stretching, const string& stretchingAccu, const string& bending);

    virtual void vectorize() override {}
    virtual void cleanup() override {}

    virtual string getName() {return "Membrane";}
    virtual void whoIsCulprit();
    
    virtual floatingpoint computeEnergy(floatingpoint* coord, bool stretched) override;
    virtual void computeForces(floatingpoint* coord, floatingpoint* f) override;
    
    virtual void computeLoadForces() { return; }
    
    virtual vector<NeighborList*> getNeighborLists() { return vector<NeighborList*>{}; }
};

#endif
