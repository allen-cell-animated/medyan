
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

#ifndef M3SYM_CController_h
#define M3SYM_CController_h

#include "common.h"

#include "ChemSim.h"
#include "ChemManager.h"

//FORWARD DECLARATIONS
class SubSystem;

/// Used to intialize, control, and run the chemical components of a simulation

/*!
 *  ChemController is a class used by the SubSystem to instantiate, control, and run
 *  the chemical dynamics of a simulation. It has functions to initialize a chemical 
 *  system, which, based on a choice of the reaction-diffusion algorithm as well as the 
 *  type of manager which controls the reactions in the simulation, as well as run the 
 *  chemical components of the simulation.
 *
 *  The controller initializes all chemical singletons used, including ChemSim
 *  and ChemManager to the correct implementations, given that they are implemented.
 */
class CController {
   
private:
    SubSystem* _subSystem; ///< A pointer to the subsystem
    
public:
    /// Constructor which sets subsystem pointer
    CController(SubSystem* s) {_subSystem = s;}
    
    /// Initialize the ChemSim algorithm as well as the ChemManager
    ///@param chemAlgorithm - a string defining the chemical algorithm to be used
    ///@param chemInitializer - a string defining the chemical manager used
    void initialize(string& chemAlgorithm, ChemistryData& chem);
    
    ///Run a number of chemical steps
    bool run(int steps) {
        
        //update copy numbers
        ChemManager::updateCopyNumbers();
        
        //run the steps
        return ChemSim::run(steps);
    }
};


#endif
