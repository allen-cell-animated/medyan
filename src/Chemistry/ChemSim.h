
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_ChemSim_h
#define MEDYAN_ChemSim_h

#include <memory>

#include "common.h"

    
//FORWARD DECLARATIONS
class ChemSimImpl;
class ReactionBase;
class DissipationTracker;

/// Used to manage running a network of chemical reactions.

/*! ChemSim implements a Strategy pattern, allowing custom algorithms for running 
 *  stochastic chemical simulations. It itself has a pointer to a  single static 
 *  implementation of ChemSimImpl. After the specific algorithm is chosen and ChemSim 
 *  is instantiated, ChemSim can be used to manage simulations, through such methods 
 *  as run(steps) etc.
 */
class ChemSim {
public:
    /// SetInstance
    /// @param ChemSimImpl *csi is a pointer the concrete implementation
    /// of the stochastic simulation algorithm.
    /// @note ChemSim simply stores the csi pointer but does not manage its memory.
    /// Make sure that csi is always a valid pointer while ChemSim is used.
    void setInstance(ChemSimImpl *csi);
    
    /// After all initial reactions have been added via addReaction(...) method,
    /// invoke initialize() prior to invoking run()
    void initialize();
    
    /// Add Reaction *r to the chemical network which needs to be simulated
    void addReaction(ReactionBase *r);
    
    /// Remove Reaction *r from the simulated chemical network 
    void removeReaction(ReactionBase *r);
    
    /// Run the chemical dynamics for a set amount of time
    bool run(double time);
    
    /// Run the chemical dynamics for a set amount of reaction steps
    bool runSteps(int steps);
    
    /// Mainly used for debugging: print chemical reactions in the network at
    /// this moment
    void printReactions();
    
    vector<double> getEnergy();
    
    DissipationTracker* getDT();
    
    
    
    
private:
    ChemSimImpl* _pimpl; ///< Store a pointer to a specific implementation
                         ///< of stochastic chemical kinetics; no ownership
};
#endif
