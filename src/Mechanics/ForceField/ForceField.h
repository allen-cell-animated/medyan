
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

#ifndef MEDYAN_ForceField_h
#define MEDYAN_ForceField_h

#include <vector>

#include "common.h"

//FORWARD DECLARATIONS
class NeighborList;
class HybridNeighborList;

/// An abstract class to represent various force field calculations
/*!
 *  ForceField is used for force calculations between elements in the SubSystem.
 *  Specific implementations of the ForceField class will have different potentials.
 */
class ForceField {
    
public:
    /// Get the name of this forcefield
    virtual string getName() = 0;
    
    /// Produce a vectorized version of interaction data
    /// Could include constants, positions, neighbors list data, etc
    virtual void vectorize() = 0;
    /// Cleanup all vectorized data
    virtual void cleanup() = 0;
    
    /// Compute total energy of this forcefield in the system
    /// @return  the energy value if valid. If an inf or NaN value has been
    /// calculated, return -1.
    virtual floatingpoint computeEnergy(floatingpoint *coord, bool stretched = false) = 0;

    /// Compute forces of this forcefield in the system. Update Bead
    /// forces accordingly.
    virtual void computeForces(floatingpoint *coord, floatingpoint *f) = 0;
    
    ///Compute all load forces on beads in this system.
    ///Updates all Bead's load force components for Reaction updating.
    virtual void computeLoadForces() = 0;
    
    /// In the case of a calculation error, print the culprit of the FF error.
    /// Typically, will just print the Trackable element where the error came from.
    virtual void whoIsCulprit() = 0;
    
    /// Get all neighbor lists associated with a ForceField
    virtual vector<NeighborList*> getNeighborLists() = 0;

    virtual void setHNeighborLists(HybridNeighborList* Hnl){};
    // assign stretchforces for Linker and Motor. Can be extended to other FFs as well.
    virtual void assignforcemags(){};
};

#endif
