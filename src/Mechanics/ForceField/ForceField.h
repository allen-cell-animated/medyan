
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

#ifndef MEDYAN_ForceField_h
#define MEDYAN_ForceField_h

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
    virtual double computeEnergy(double d, double *coord) = 0;
    /// Compute forces of this forcefield in the system. Update Bead
    /// forces accordingly.
    virtual void computeForces(double *coord, double *f) = 0;
    /// Compute auxiliary forces of this forcefield in the system. Update
    /// Bead auxiliary forces accordingly.
    virtual void computeForcesAux(double *coord, double *fa) = 0;
    
    ///Compute all load forces on beads in this system.
    ///Updates all Bead's load force components for Reaction updating.
    virtual void computeLoadForces() = 0;
    
    /// In the case of a calculation error, print the culprit of the FF error.
    /// Typically, will just print the Trackable element where the error came from.
    virtual void whoIsCulprit() = 0;
    
    /// Get all neighbor lists associated with a ForceField
    virtual vector<NeighborList*> getNeighborLists() = 0;
    
    
};

#endif
