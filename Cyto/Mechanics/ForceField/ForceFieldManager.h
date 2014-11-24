//
//  ForceFieldManager.h
//  Cyto
//
//  Created by James Komianos on 9/10/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__ForceFieldManager__
#define __Cyto__ForceFieldManager__

#include <iostream>
#include <vector>

#include "common.h"

#include "ForceField.h"
#include "BeadDB.h"

///ForceFieldManager is a class to store and iterate over all forcefields.
/*!
 *  The ForceFieldManager is used to store all forcefields initialized by the
 *  system, as well as iterate over these forcefields and calculate total
 *  forces and energies. Contains functions for the said calculations.
 */
class ForceFieldManager {
    
public:
     vector<ForceField*> _forceFields;
    
    //Compute the energy using all available force fields
    double computeEnergy(double d) {
        
        double energy = 0;
        for(auto &f : _forceFields)
            energy += f->computeEnergy(d);
        /// pass it to subsystem!!!
        return energy;
    }
    
    ///Compute the forces of all force fields
    void computeForces() {
        resetForces();
        for(auto &f : _forceFields) f->computeForces();
    }
    
    ///Compute the forcesAux of all force fields
    void computeForcesAux() {
        resetForcesAux();
        
        for(auto &f : _forceFields)
            f->computeForcesAux();
    }
    
    ///Reset the forces of all objects
    void resetForces() {
        
        for(auto it: *BeadDB::instance()) {
            it->force.assign (3, 0); //Set force to zero;
        }
    }
    
    ///Reset the forcesAux of all objects
    void resetForcesAux() {
        
        for(auto it: *BeadDB::instance()) {
            it->forceAux.assign (3, 0); //Set forceAux to zero;
        }
    }
    
};



#endif /* defined(__Cyto__ForceFieldManager__) */
