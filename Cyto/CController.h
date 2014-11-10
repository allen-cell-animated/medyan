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
#include "common.h"
#include "SubSystem.h"
#include "ChemSim.h"
#include "ChemManager.h"

#include "ChemNRMImpl.h"
#include "ChemGillespieImpl.h"
#include "ChemSimpleGillespieImpl.h"
#include "SimpleManagerImpl.h"

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
    void initialize(string& chemAlgorithm, string chemInitializer, ChemistryData& chem);
    
    ///Run a number of chemical steps
    bool run(int steps) { return ChemSim::run(ChemSimRunKey(), steps); }
};
    


#endif /* defined(__Cyto__CController__) */
