
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

#ifndef MEDYAN_CamkiiFF_h
#define MEDYAN_CamkiiFF_h

#include <vector>

#include "common.h"

#include "ForceField.h"

//FORWARD DECLARATIONS
class CamkiiInteractions;
class Camkii;

/// An implementation of the ForceField class that calculates Filament interactions.
class CamkiiFF : public ForceField {
 
private:
    vector<unique_ptr<CamkiiInteractions>>
    _interactions; ///< Vector of initialized filament interactions
    
    CamkiiInteractions* _culpritInteraction; ///< Culprit in case of error
public:
    /// Constructor, intializes stretching, bending, and twisting forces
    CamkiiFF(string& stretching, string& bending);
    
    virtual string getName() {return "Camkii";}
    virtual void whoIsCulprit();
    
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();
    
    virtual void computeLoadForces() {return;}
    
    virtual vector<NeighborList*> getNeighborLists() {return vector<NeighborList*>{};}
};

#endif
