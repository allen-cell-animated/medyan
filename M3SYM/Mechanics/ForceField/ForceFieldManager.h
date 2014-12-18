
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

#ifndef M3SYM_ForceFieldManager_h
#define M3SYM_ForceFieldManager_h

#include <vector>

#include "common.h"

#include "ForceField.h"
#include "Bead.h"

/// A class to store and iterate over all [ForceFields](@ref ForceField).
/*!
 *  The ForceFieldManager is used to store all [ForceFields](@ref ForceField) initialized by the
 *  system, as well as iterate over these potentials and calculate total
 *  forces and energies. This class contains functions for the said calculations.
 */
class ForceFieldManager {
    
public:
     vector<ForceField*> _forceFields; ///< All forcefields in the system
    
    /// Compute the energy using all available force fields
    double computeEnergy(double d) {
        
        double energy = 0;
        for(auto &f : _forceFields)
            energy += f->computeEnergy(d);
        
        //if energy is infinity, exit ungracefully.
        if(energy == numeric_limits<double>::infinity()) {
            cout <<
            "Energy became infinite. Try adjusting equilibration step size." << endl;
            exit(EXIT_FAILURE);
        }
        return energy;
    }
    
    /// Compute the forces of all force fields
    void computeForces() {
        resetForces();
        for(auto &f : _forceFields) f->computeForces();
    }
    
    /// Compute the forcesAux of all force fields
    void computeForcesAux() {
        resetForcesAux();
        for(auto &f : _forceFields)
            f->computeForcesAux();
    }
    
    /// Reset the forces of all objects
    void resetForces() {
        
        for(auto it: *BeadDB::instance())
            it->force.assign (3, 0); //Set force to zero;
    }
    
    /// Reset the forcesAux of all objects
    void resetForcesAux() {
        
        for(auto it: *BeadDB::instance())
            it->forceAux.assign (3, 0); //Set forceAux to zero;
    }
};

#endif
