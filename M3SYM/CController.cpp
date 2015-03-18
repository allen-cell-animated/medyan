
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "CController.h"

#include "SubSystem.h"

#include "ChemNRMImpl.h"
#include "ChemGillespieImpl.h"
#include "ChemSimpleGillespieImpl.h"
#include "SimpleManagerImpl.h"

void CController::initialize(string& chemAlgorithm, string& chemManager,
                             ChemistryData& chem) {
    
    // Set instance of chemsim algorithm
    ChemSimImpl* csi;
    if(chemAlgorithm == "NRM") {
        
#if !defined(TRACK_DEPENDENTS)
        cout << "The NRM algorithm relies on tracking dependents. Please set this"
            << " compilation macro and try again. Exiting." << endl;
        exit(EXIT_FAILURE);
#endif
        csi = new ChemNRMImpl;
    }
    
    else if(chemAlgorithm == "GILLESPIE") {
        
#if !defined(TRACK_DEPENDENTS)
        cout << "The Gillespie algorithm relies on tracking dependents. Please set this"
            << " compilation macro and try again. Exiting." << endl;
        exit(EXIT_FAILURE);
#endif
        csi = new ChemGillespieImpl;
    }
    
    else if(chemAlgorithm == "SIMPLEGILLESPIE")
        csi = new ChemSimpleGillespieImpl;
    
    else {
        cout<< "Chem algorithm not recognized. Exiting." <<endl;
        exit(EXIT_FAILURE);
    }
    ChemSim::setInstance(csi);
    
    // Set the instance of the initializer
    ChemManagerImpl* cii;
    cii = new SimpleManagerImpl(_subSystem);

    ChemManager::setInstance(cii);
    ChemManager::initialize(chem);
    ChemSim::initialize();
}

