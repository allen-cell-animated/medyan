
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_CController_h
#define MEDYAN_CController_h

#include "common.h"

#include "Compartment.h"

//FORWARD DECLARATIONS
class SubSystem;
class ChemSim;
class ChemManager;
class CompartmentGrid;

struct ChemistryData;

/// Used to intialize, control, and run the chemical components of a simulation

/*!
 *  ChemController is a class used by the SubSystem to instantiate, control, and run
 *  the chemical dynamics of a simulation. It has functions to initialize a chemical 
 *  system, which, based on a choice of the reaction-diffusion algorithm as well as the 
 *  type of manager which controls the reactions in the simulation, as well as run the 
 *  chemical components of the simulation.
 *
 *  The controller initializes all chemical objects used, including ChemSim
 *  and ChemManager to the correct implementations, given that they are implemented.
 */
class CController {
   
private:
    SubSystem* _subSystem; ///< A pointer to the SubSystem
    
    //@{
    /// Holding pointers to control all chemistry.
    ChemSim* _chemSim; ChemManager* _chemManager;
    //@}
    
public:
    /// Constructor which sets subsystem pointer
    CController(SubSystem* s) : _subSystem(s) {}
    
    ///Activate a compartment. Wrapper function for Compartment::activate().
    void activate(Compartment* C, Compartment::ActivateReason reason = Compartment::ActivateReason::Whole) { C->activate(_chemSim, reason); }
    /// Update activation of a compartment. Wrapper function for Compartment::updateActivation()
    void updateActivation(Compartment* C, Compartment::ActivateReason reason = Compartment::ActivateReason::Whole) { C->updateActivation(_chemSim, reason); }
    ///Deactivate a compartment. Wrapper function for Compartment::deactivate().
    void deactivate(Compartment* C, bool init=false) {C->deactivate(_chemSim, init);}
    
    /// Initialize the ChemSim algorithm as well as the ChemManager
    ///@param chemAlgorithm - a string defining the chemical algorithm to be used
    ///@param chemInitializer - a string defining the chemical manager used
    void initialize(string& chemAlgorithm, ChemistryData& chem);
    //aravind June 29,2016.
    void restart();
    
    ///Run chemistry for a given amount of time
    bool run(double time);
    
    ///Run chemistry for a given number of reaction steps
    bool runSteps(int steps);
    
    ///Remove set of reactions at runtime, specified by input
    void removeReactions();
};


#endif
