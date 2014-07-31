//
//  ChemController.h
//  Cyto
//
//  Created by James Komianos on 7/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__ChemController__
#define __Cyto__ChemController__

#include <iostream>
#include "ChemSim.h"
#include "ChemInitializer.h"

#include "ChemNRMImpl.h"
#include "ChemGillespieImpl.h"
#include "ChemSimpleGillespieImpl.h"
#include "SimpleInitializerImpl.h"

namespace chem {


    ///ChemController class is used to intialize, control, and run the chemical components of a simulation
    
    /*!
     *  ChemController is a class used by the SubSystem class to instantiate, control, and run the chemical
     *  dynamics of a simulation. It has functions to initialize a chemical system, which, based on a choice
     *  of the reaction-diffusion algorithm as well as the type of initializer which controls the reactions
     *  in the simulation, as well as run the chemical components of the simulation.
     *
     *  The controller initializes all chemical singletons used, including ChemSim and ChemInitializer,
     *  to the correct implementations, given that they are implemented.
     */
    class ChemController {
       
        ///Initialize the chemical system. MUST BE CALLED BEFORE RUN!
        
        ///@param chemAlgorithm - a string defining the chemical algorithm to be used
        ///       could be: NRM, Gillespie, SimpleGillespie.
        ///@param chemInitializer - a string defining the chemical initializer used
        ///       could be: Simple.
        void initialize(std::string chemAlgorithm, std::string chemInitializer) {
            
            ///Set instance of chemsim algorithm
            ChemSimImpl* csi;
            if(chemAlgorithm == "NRM") {
                ChemNRMImpl chem;
                csi = &chem;
            }
            else if(chemAlgorithm == "Gillespie") {
                ChemGillespieImpl chem;
                csi = &chem;
            }
            else if(chemAlgorithm == "SimpleGillespie") {
                ChemSimpleGillespieImpl chem;
                csi = &chem;
            }
            else {
                std::cout<< "Chem algorithm not found. Exiting." <<std::endl;
                exit(EXIT_FAILURE);
            }
            ChemSim::setInstance(ChemSimInitKey(), csi);
            
            
            ///Set the instance of the initializer
            ChemInitializerImpl* cii;
            if(chemInitializer == "Simple") {
                SimpleInitializerImpl init;
                cii = &init;
            }
            else {
                std::cout<< "Initializer type not found. Exiting." <<std::endl;
                exit(EXIT_FAILURE);
            }
            ChemInitializer::setInstance(ChemInitializerInitKey(), cii);
            
            
            
            
            ///initialize grid ...
            //ChemInitializer::initializeGrid(g);
            
            
            
            
            ///initialize chemsim
            ChemSim::initialize(ChemSimInitKey());
            
            
            
        }
    };

        
        
} ///end namespace chem
        
        
        
        
    



















#endif /* defined(__Cyto__ChemController__) */
