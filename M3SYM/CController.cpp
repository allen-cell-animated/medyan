
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

#include "ChemManager.h"

#include "ChemNRMImpl.h"
#include "ChemGillespieImpl.h"
#include "ChemSimpleGillespieImpl.h"

#include "CCylinder.h"
#include "Cylinder.h"

void CController::initialize(string& chemAlgorithm, ChemistryData& chem) {
    
    // new ChemSim object
    _chemSim = new ChemSim;
    
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
    _chemSim->setInstance(csi);
    
    //Create manager, intialize
    _chemManager = new ChemManager(_subSystem, chem);
    _chemManager->initializeSystem(_chemSim);
    
    // init chemsim
    _chemSim->initialize();
    
    // set some static ptrs
    FilamentReactionTemplate::_ps = _subSystem;
    
    CCylinder::_chemSim = _chemSim;
    Cylinder::_chemManager = _chemManager;
    
}

bool CController::run(double time) {
    
    //update copy numbers
    _chemManager->updateCopyNumbers();
    
    //run the steps
    return _chemSim->run(time);
}

