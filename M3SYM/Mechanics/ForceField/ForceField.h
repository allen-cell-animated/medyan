
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

#ifndef M3SYM_ForceField_h
#define M3SYM_ForceField_h

#include "common.h"

/// An abstract class to represent various force field calculations
/*!
 *  ForceField is used for force calculations between filaments, beads, linkers, etc.
 *  Specific implementations of the ForceField class will have different potentials.
 */
class ForceField {

private:
    string _name; ///< Name of this forcefield
    
public:
    /// Get the name of this forcefield
    const string& getName() {return _name;}
    /// Compute total energy of this forcefield in the system
    virtual double computeEnergy(double d) = 0;
    /// Compute forces of this forcefield in the system. Update [Bead](@ref Bead) forces accordingly.
    virtual void computeForces() = 0;
    /// Compute auxiliary forces of this forcefield in the system. Update [Bead](@ref Bead) auxiliary forces accordingly.
    virtual void computeForcesAux() = 0;
};

#endif
