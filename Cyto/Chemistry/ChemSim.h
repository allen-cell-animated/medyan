
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_ChemSim_h
#define M3SYM_ChemSim_h

#include <memory>

#include "common.h"
#include "utility.h"
    
//FORWARD DECLARATIONS
class ChemSimImpl;
class ReactionBase;

/// ChemSim is used to manage running a network of chemical reactions.

/*! ChemSim implements a Strategy pattern, allowing custom algorithms for running stochastic chemical simulations. 
 *  It itself is a singleton, with a pointer to a single static implementation of ChemSimImpl.
 *  After the specific algorithm is chosen and ChemSim is instantiated, ChemSim can be used to manage simulations, through 
 *  such methods as run(steps) etc. Here is an example (assuming key is valid for all functions):
 *  @code
        SpeciesBulk A1("A1",  25);
        SpeciesBulk A2("A2", 25);
        Reaction r1 = { {&A1,&A2}, 1, 1, 100.0 };
        ChemSim::setInstance(key, new ChemNRMImpl());
        ChemSim::addReaction(key, &r1);
        ChemSim::initialize(key);
        ChemSim::run(key)
 * @endcode
 *
 *  Specific functions in this class require a key to be accessed. The key declarations are above. These keys can
 *  only be created and/or destroyed by the classes that are "friends" with the key.
 */
class ChemSim {
public:
    /// SetInstance
    /// @param ChemSimImpl *csi is a pointer the concrete implementation of stochastic simulation algorithm.
    /// @note ChemSim simply stores the csi pointer but does not manage its memory. Make sure that csi is always 
    /// a valid pointer while ChemSim is used.
    static void setInstance(ChemSimImpl *csi);
    
    /// Copying is not allowed
    ChemSim(const ChemSim &rhs) = delete;

    /// Assignment is not allowed
    ChemSim& operator=(ChemSim &rhs) = delete;
    
    /// After all initial reactions have been added via addReaction(...) method, invoke initialize() prior to invoking run() 
    static void initialize();
    
    /// Add Reaction *r to the chemical network which needs to be simulated
    static void addReaction(ReactionBase *r);
    
    /// Remove Reaction *r from the simulated chemical network 
    static void removeReaction(ReactionBase *r);
    
    /// Run the chemical dynamics for a specific number of steps
    static bool run(int steps);
    
    /// Mainly used for debugging: print chemical reactions in the network at this moment
    static void printReactions();
private:

    static ChemSimImpl* _pimpl; ///< Store a pointer to a specific implementation of stochastic chemical kinetics; no ownership
    ChemSim() {};
};
#endif
