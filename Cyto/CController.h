//
//  CController.h
//  Cyto
//
//  Created by James Komianos on 7/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__CController__
#define __Cyto__CController__

#include <iostream>
#include "SubSystem.h"
#include "ChemSim.h"
#include "ChemInitializer.h"

#include "ChemNRMImpl.h"
#include "ChemGillespieImpl.h"
#include "ChemSimpleGillespieImpl.h"
#include "SimpleInitializerImpl.h"

///CController class is used to intialize, control, and run the chemical components of a simulation

/*!
 *  ChemController is a class used by the SubSystem class to instantiate, control, and run the chemical
 *  dynamics of a simulation. It has functions to initialize a chemical system, which, based on a choice
 *  of the reaction-diffusion algorithm as well as the type of initializer which controls the reactions
 *  in the simulation, as well as run the chemical components of the simulation.
 *
 *  The controller initializes all chemical singletons used, including ChemSim and ChemInitializer,
 *  to the correct implementations, given that they are implemented.
 */
class CController {
   
private:
    SubSystem* _subSystem;
    
public:
    CController(SubSystem* s) {_subSystem = s;}
    ///Initialize the chemical system. MUST BE CALLED BEFORE RUN!
    
    ///@param chemAlgorithm - a string defining the chemical algorithm to be used
    ///       could be: NRM, Gillespie, SimpleGillespie.
    ///@param chemInitializer - a string defining the chemical initializer used
    ///       could be: Simple.
    void initialize(std::string& chemAlgorithm, std::string& chemInitializer) {
        
        ///Set instance of chemsim algorithm
        ChemSimImpl* csi;
        if(chemAlgorithm == "NRM")
            csi = new ChemNRMImpl;
        
        else if(chemAlgorithm == "GILLESPIE")
            csi = new ChemGillespieImpl;
        
        else if(chemAlgorithm == "SIMPLEGILLESPIE")
            csi = new ChemSimpleGillespieImpl;
        
        else {
            std::cout<< "Chem algorithm not found. Exiting." <<std::endl;
            exit(EXIT_FAILURE);
        }
        ChemSim::setInstance(ChemSimInitKey(), csi);
        
        
        ///Set the instance of the initializer
        ChemInitializerImpl* cii;
        if(chemInitializer == "SIMPLE") {
            cii = new SimpleInitializerImpl;
        }
        else {
            std::cout<< "Initializer type not found. Exiting." <<std::endl;
            exit(EXIT_FAILURE);
        }
        ChemInitializer::setInstance(ChemInitializerInitKey(), cii);
        
        ///initialize grid ...
        ChemInitializer::initializeGrid(ChemInitializerGridKey());

        ///initialize chemsim
        ChemSim::initialize(ChemSimInitKey());
    }
    
    ///Run a number of chemical steps
    bool run(int steps) { return ChemSim::run(ChemSimRunKey(), steps); }
};
    


#endif /* defined(__Cyto__CController__) */