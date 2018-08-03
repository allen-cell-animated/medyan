
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

#ifndef MEDYAN_CaMKIIingFF_h
#define MEDYAN_CaMKIIingFF_h

#include <vector>

#include "common.h"

#include "ForceField.h"

//FORWARD DECLARATIONS
class CaMKIIingInteractions;
class CaMKIIingPoint;

/// CaMKIIing FF is an implementation of the ForceField class that
/// calculates CaMKIIingPoint interactions.
class CaMKIIingFF : public ForceField {
    
private:
    vector<unique_ptr<CaMKIIingInteractions>>
    _camkiiingInteractionVector; ///< Vector of initialized camkiiing interactions
    
    CaMKIIingInteractions* _culpritInteraction; ///< Culprit in case of error
public:
    /// Constructor, intializes all interaction at the camkiiing point
    CaMKIIingFF(string& stretching, string& bending,
                string& dihedral, string& position);
    
    virtual string getName() {return "CaMKIIing";}
    virtual void whoIsCulprit();
    
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();
    
    virtual void computeLoadForces() {return;}
    
    virtual vector<NeighborList*> getNeighborLists() {return vector<NeighborList*>{};}
};

#endif
