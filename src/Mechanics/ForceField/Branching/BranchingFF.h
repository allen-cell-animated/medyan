
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_BranchingFF_h
#define MEDYAN_BranchingFF_h

#include <vector>

#include "common.h"

#include "ForceField.h"

//FORWARD DECLARATIONS
class BranchingInteractions;
class BranchingPoint;

/// Branching FF is an implementation of the ForceField class that
/// calculates BranchingPoint interactions.
class BranchingFF : public ForceField {
    
private:
    vector<unique_ptr<BranchingInteractions>>
    _branchingInteractionVector; ///< Vector of initialized branching interactions
    
protected:
    BranchingInteractions* _culpritInteraction; ///< Culprit in case of error
    
public:
    /// Constructor, intializes all interaction at the branching point
    BranchingFF(string& stretching, string& bending,
                string& dihedral, string& position);
    
    virtual void vectorize();
    virtual void cleanup();
    
    virtual string getName() {return "Branching";}
    virtual void whoIsCulprit();
    
    virtual double computeEnergy(double *coord, double *f, double d);
    virtual void computeForces(double *coord, double *f);
    
    virtual void computeLoadForces() {return;}
    
    virtual vector<NeighborList*> getNeighborLists() {return vector<NeighborList*>{};}
};

#endif
