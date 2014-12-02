
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

#ifndef M3SYM_MotorGhostFF_h
#define M3SYM_MotorGhostFF_h

#include <vector>

#include "common.h"

#include "ForceField.h"

//FORWARD DECLARATIONS
class MotorGhostInteractions;

/// An implementation of the [ForceField](@ref ForceField) class that calculates [Motor] (@ref Motor)
/// stretching, bending, and twisting.
class MotorGhostFF : public ForceField {
    
private:
    vector <unique_ptr<MotorGhostInteractions>> _motorGhostInteractionVector; ///< Vector of initialized motor interactions
    
public:
    /// Constructor, intializes stretching, bending, and twisting forces
    MotorGhostFF(string& stretching, string& bending, string& twisting);
    
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();
};

#endif
