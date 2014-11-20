//
//  CController.cpp
//  Cyto
//
//  Created by James Komianos on 7/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "CController.h"

void CController::initialize(string& chemAlgorithm, string chemInitializer, ChemistryData& chem) {
    
    ///Set instance of chemsim algorithm
    ChemSimImpl* csi;
    if(chemAlgorithm == "NRM")
        csi = new ChemNRMImpl;
    
    else if(chemAlgorithm == "GILLESPIE")
        csi = new ChemGillespieImpl;
    
    else if(chemAlgorithm == "SIMPLEGILLESPIE")
        csi = new ChemSimpleGillespieImpl;
    
    else {
        cout<< "Chem algorithm not found. Exiting." <<endl;
        exit(EXIT_FAILURE);
    }
    ChemSim::setInstance(csi);
    
    
    ///Set the instance of the initializer
    ChemManagerImpl* cii;
    cii = new SimpleManagerImpl(_subSystem);


    ChemManager::setInstance(cii);
    
    ///initialize grid ...
    ChemManager::initialize(chem);
    
    ///initialize chemsim
    ChemSim::initialize();
}

