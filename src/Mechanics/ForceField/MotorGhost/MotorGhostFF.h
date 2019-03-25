
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

#ifndef MEDYAN_MotorGhostFF_h
#define MEDYAN_MotorGhostFF_h

#include <vector>

#include "common.h"

#include "ForceField.h"

//FORWARD DECLARATIONS
class MotorGhostInteractions;
class MotorGhost;

/// An implementation of the ForceField class that calculates MotorGhost
/// stretching, bending, and twisting.
class MotorGhostFF : public ForceField {
    
private:
    vector <unique_ptr<MotorGhostInteractions>>
    _motorGhostInteractionVector; ///< Vector of initialized motor interactions
    
protected:
    MotorGhostInteractions* _culpritInteraction; ///< Culprit in case of error
    
public:
    /// Constructor, intializes stretching, bending, and twisting forces
    MotorGhostFF(string& stretching, string& bending, string& twisting);
    
    virtual void vectorize();
    virtual void cleanup();

    virtual string getName() {return "MotorGhost";}
    virtual void whoIsCulprit();
    
    virtual totalenergyfloatingpoint computeEnergy(floatingpoint *coord, totalforcefloatingpoint *f, floatingpoint d);
    virtual void computeForces(floatingpoint *coord, totalforcefloatingpoint *f);
    
    virtual void computeLoadForces() {return;}
    
    virtual vector<NeighborList*> getNeighborLists() {return vector<NeighborList*>{};}
    //Assigns stretchforces for ratechangeimpl
    virtual void assignforcemags();
};

#endif
