
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
    template< bool stretched > void updateGeometryValue() const;
    void updateGeometryValueWithDerivative() const;

    /// Compute the energy using all available force fields
    /// @return Returns infinity if there was a problem with a ForceField
    /// energy calculation, such that beads will not be moved to this
    /// problematic configuration.
    /// @param print - prints detailed info about energies
    template< bool stretched = false > double computeEnergy(bool verbose = false) {
        double energy = 0;
        for(auto &f : _forceFields) {
            
            auto tempEnergy = f->computeEnergy(stretched);
            
            if(verbose) cout << f->getName() << " energy = " << tempEnergy << endl;
            
            //if energy is infinity, exit with infinity.
            if(tempEnergy <= -1) {
                
                //if this is the current energy, exit ungracefully
                if(!stretched) {
                    
                    cout << "Energy = " << tempEnergy << endl;
                    
                    cout << "Energy of system became infinite. Try adjusting minimization parameters." << endl;
                    cout << "The culprit was ... " << f->getName() << endl;
                    
                    //get the culprit in output
                    f->whoIsCulprit();
                    
                    exit(EXIT_FAILURE);
                }
                //if this is a minimization try, just return infinity
                else return numeric_limits<double>::infinity();
            }
            else energy += tempEnergy;
        }
        return energy;

    }
    
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
