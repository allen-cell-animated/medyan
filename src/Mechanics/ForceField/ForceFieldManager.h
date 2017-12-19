
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

#ifndef MEDYAN_ForceFieldManager_h
#define MEDYAN_ForceFieldManager_h

#include <vector>

#include "common.h"

#include "ForceField.h"
#include "Bead.h"

// Forward declarations
class SubSystem;

/// A class to store and iterate over all [ForceFields](@ref ForceField).
/*!
 *  The ForceFieldManager is used to store all [ForceFields](@ref ForceField) 
 *  initialized by the system, as well as iterate over these potentials and calculate 
 *  total forces and energies. This class contains functions for the said calculations.
 */
class ForceFieldManager {
    
public:
    vector<ForceField*> _forceFields; ///< All forcefields in the system

    SubSystem* _subSystem; ///< Pointer to the subsystem

    ForceFieldManager(SubSystem* s) { _subSystem = s; }
    
    /// Update the geometry of all elements in the system
    void updateGeometries(bool calcDerivative=false, double d=0.0);
    
    /// Compute the energy using all available force fields
    /// @return Returns infinity if there was a problem with a ForceField
    /// energy calculation, such that beads will not be moved to this
    /// problematic configuration.
    /// @param print - prints detailed info about energies
    double computeEnergy(double d, bool verbose = false);
    
    /// Compute the forces of all force fields 
    void computeForces();
    /// Compute the forcesAux of all force fields
    void computeForcesAux();
    /// Compute forcesAuxP of all force fields
    void computeForcesAuxP();
    
    /// Compute the load forces on the beads. This does not update the force (xyz) vector
    /// contained by Bead, but updates the loadForce vector which contains precalculated
    /// load values based on the bead's directionality of growth in a filament.
    void computeLoadForces();
    
    /// Reset the forces of all objects
    void resetForces();
    /// Reset the forcesAux of all objects
    void resetForcesAux();
};

#endif
