//
//  MController.h
//  Cyto
//
//  Created by James Komianos on 8/4/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__MController__
#define __Cyto__MController__

#include "SubSystem.h"
#include "ForceField.h"
#include "Mcommon.h"
#include <iostream>

/// MController class is used to initialize and run the mechanical components of a simulation

/*!
 *  MechanicalController is a class used by the SubSystem to initialize force fields, given an initial
 *  selection of which force fields should be included. It can compute forces and energies
 *  over all force fields, as well as run energy minimization algorithms.
 */

class MController {
    
private:
    std::vector<ForceField> _forceFields; ///< vector of force field selections
    
public:
    
    ///Initialize the MController using a list of vector names
    ///@param forceFields - a list of forcefields to be added
    void initialize (std::initializer_list<std::string> forceFields)
    {
        
        for(auto &f : forceFields) {
            ///implement this
        }
    }
    
    ///Compute the energy using all available force fields
    double ComputeEnergy() {
        
        double energy = 0;
        for(auto &f : _forceFields)
            energy += f.ComputeEnergy();
        
        return energy;
    }
    
    ///Reset the forces of all objects
    void ResetForces() {
        
        ///implement this
    }
    
    ///Compute the forces of all force fields
    void ComputeForces() {
        ResetForces();
        
        for(auto &f : _forceFields)
            f.ComputeForces();
    }
    
    
    ///Run minimization on the system using the chosen algorithm
   void run(System* ps, std::string solver) {
     
        if(solver == "FletcherRieves")
            FletcherRievesMethod(ps);
        else if (solver == "PolakRibiere")
            PolakRibiereMethod(ps);
        
        else {
            std::cout<< "Mechanical algorithm not found. Exiting." <<std::endl;
            exit(EXIT_FAILURE);
        }
    }
    
    
    
};




#endif /* defined(__Cyto__MController__) */
