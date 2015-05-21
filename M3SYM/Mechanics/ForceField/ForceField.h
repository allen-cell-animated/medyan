
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

#ifndef M3SYM_ForceField_h
#define M3SYM_ForceField_h

#include <vector>

#include "common.h"

//FORWARD DECLARATIONS
class NeighborList;

/// An abstract class to represent various force field calculations
/*!
 *  ForceField is used for force calculations between elements in the SubSystem.
 *  Specific implementations of the ForceField class will have different potentials.
 */
class ForceField {
    
public:
    /// Get the name of this forcefield
    virtual string getName() = 0;
    
    /// Compute total energy of this forcefield in the system
    /// @return  the energy value if valid. If an inf or NaN value has been
    /// calculated, return -1.
    virtual double computeEnergy(double d) = 0;
    /// Compute forces of this forcefield in the system. Update Bead
    /// forces accordingly.
    virtual void computeForces() = 0;
    /// Compute auxiliary forces of this forcefield in the system. Update
    /// Bead auxiliary forces accordingly.
    virtual void computeForcesAux() = 0;
    
    /// Get all neighbor lists associated with a ForceField
    virtual vector<NeighborList*> getNeighborLists() = 0;
};

#endif
