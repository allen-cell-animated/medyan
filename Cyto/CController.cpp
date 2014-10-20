//
//  CController.cpp
//  Cyto
//
//  Created by James Komianos on 7/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "CController.h"

void CController::initialize(std::string& chemAlgorithm, std::string chemInitializer, ChemistrySpeciesAndReactions& chemSR) {
    
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
        cii = new SimpleInitializerImpl(_subSystem);
    }
    else {
        std::cout<< "Initializer type not found. Exiting." <<std::endl;
        exit(EXIT_FAILURE);
    }
    ChemInitializer::setInstance(ChemInitializerInitKey(), cii);
    
    ///initialize grid ...
    ChemInitializer::initialize(ChemInitializerGridKey(), chemSR);
    
    ///initialize chemsim
    ChemSim::initialize(ChemSimInitKey());
}

