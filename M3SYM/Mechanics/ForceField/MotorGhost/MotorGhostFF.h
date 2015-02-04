
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_MotorGhostFF_h
#define M3SYM_MotorGhostFF_h

#include <vector>

#include "common.h"

#include "ForceField.h"

//FORWARD DECLARATIONS
class MotorGhostInteractions;

/// An implementation of the ForceField class that calculates MotorGhost
/// stretching, bending, and twisting.
class MotorGhostFF : public ForceField {
    
private:
    vector <unique_ptr<MotorGhostInteractions>>
        _motorGhostInteractionVector; ///< Vector of initialized motor interactions
    
public:
    /// Constructor, intializes stretching, bending, and twisting forces
    MotorGhostFF(string& stretching, string& bending, string& twisting);
    
    virtual string getName() {return "MotorGhost";}
    
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();
};

#endif
