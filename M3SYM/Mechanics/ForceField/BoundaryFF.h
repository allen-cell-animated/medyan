
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_BoundaryFF_h
#define M3SYM_BoundaryFF_h

#include <vector>

#include "common.h"

#include "ForceField.h"

//FORWARD DECLARATIONS
class BoundaryInteractions;

/// An implementation of the ForceField class that calculates BoundaryElement
/// repulsion and attraction to [Beads](@ref Bead) in the system.
class BoundaryFF : public ForceField {
    
private:
    vector<unique_ptr<BoundaryInteractions>> _BoundaryInteractionVector; ///< Vector of initialized boundary element interactions
    
public:
    /// Initialize the forcefields (repulsion, attraction, etc)
    BoundaryFF(string interaction1, string interaction2, string interaction3);
    
    virtual string getName() {return "Boundary";}
    
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();
};

#endif
