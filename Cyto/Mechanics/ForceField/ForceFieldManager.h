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
#include "ForceField.h"

///ForceFieldManager is a class to store and iterate over all forcefields.
/*!
 *  The ForceFieldManager is used to store all forcefields initialized by the
 *  system, as well as iterate over these forcefields and calculate total
 *  forces and energies. Contains functions for the said calculations.
 */
class ForceFieldManager {
    
public:
     std::vector<ForceField*> _forceFields;
    
    //Compute the energy using all available force fields
    double ComputeEnergy(double d) {
        
        double energy = 0;
        for(auto &f : _forceFields)
            energy += f->ComputeEnergy(d);
        /// pass it to subsystem!!!
        return energy;
    }
    
    ///Compute the forces of all force fields
    void ComputeForces() {
        ResetForces();
        
        for(auto &f : _forceFields)
            f->ComputeForces();
    }
    
    ///Compute the forcesAux of all force fields
    void ComputeForcesAux() {
        ResetForces();
        
        for(auto &f : _forceFields)
            f->ComputeForcesAux();
    }
    
    ///Reset the forces of all objects
    void ResetForces() {
        
        ///implement this
    }
    
    ///Reset the forcesAux of all objects
    void ResetForcesAux() {
        
        ///implement this
    }
    
    
};



#endif /* defined(__Cyto__ForceFieldManager__) */
