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
#include "MMinimizer.h"
#include "Parser.h"
#include "Mcommon.h"
#include "ForceField.h"
#include "MFilamentFF.h"
#include "MLinkerFF.h"
#include "MMotorGhostFF.h"
#include "MConjugateGradient.h"
#include "MCGPolakRibiereMethod.cpp"
#include "MCGFletcherRievesMethod.h"
#include <iostream>
#include <vector>


/// MController class is used to initialize and run the mechanical components of a simulation

/*!
 *  MechanicalController is a class used by the SubSystem to initialize force fields, given an initial
 *  selection of which force fields should be included. It can compute forces and energies
 *  over all force fields, as well as run energy minimization algorithms.
 */

class MController {
    
private:
    std::vector<ForceField*> _forceFields; ///< vector of force field selections
    std::vector<Minimizer*> _minimizerAlgorithms; ///<vector with algorythms for system equlibration
    
    ///Initialize the MController using a list of vector names
    ///@param forceFields - a list of forcefields to be added
    void initializeFF (MechanicsFFType forceFields)
    {
        /// Check if exist!!!
        _forceFields.push_back(new FilamentFF(forceFields.FStretchingType, forceFields.FBendingType, forceFields.FTwistingType) );
        _forceFields.push_back(new LinkerFF(forceFields.LStretchingType, forceFields.LBendingType, forceFields.LTwistingType) );
        _forceFields.push_back(new MotorGhostFF(forceFields.MStretchingType, forceFields.MBendingType, forceFields.MTwistingType) );
        
        /// Add other FF's
        
    }
    
    void initializeMinAlorythms (MechanicsAlgorithm Minimizers)
    {
        if (Minimizers.algorithm == "FLETCHERRIVES") {_minimizerAlgorithms.push_back(new ConjugateGradient<FletcherRieves>() );}
        if (Minimizers.algorithm == "POLAKRIBIERE") {_minimizerAlgorithms.push_back(new ConjugateGradient<PolakRibiere>() );}
        
        ///Change .algorithm to .cgFR , .md etc
    }

public:
    
    void initialize(MechanicsFFType forceFields, MechanicsAlgorithm Minimizers )
    {
        initializeFF(forceFields);
        initializeMinAlorythms(Minimizers);
    }
    
    //Compute the energy using all available force fields
    double ComputeEnergy(double d) {
        
        double energy = 0;
        for(auto &f : _forceFields)
            energy += f.ComputeEnergy(d);
        
        return energy;
    }
    
    ///Reset the forces of all objects
    void ResetForces() {
        
        ///implement this
    }
    
    ///Reset the forcesAux of all objects
    void ResetForcesAux() {
        
        ///implement this
    }

    
    ///Compute the forces of all force fields
    void ComputeForces() {
        ResetForces();
        
        for(auto &f : _forceFields)
            f.ComputeForces();
    }
    
    ///Compute the forcesAux of all force fields
    void ComputeForcesAux() {
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
