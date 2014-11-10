//
//  ChemSim.h
//  CytoSim
//
//  Created by Garegin Papoian on 5/2/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_ChemSim_h
#define CytoSim_ChemSim_h

#include <memory>

#include "common.h"
#include "Reaction.h"
    
///forward delcarations
class ChemSimImpl;

///Key for initialization of ChemSim
class ChemSimInitKey { friend class CController;
#ifdef TESTING
                       public:
#endif //TESTING
                       ChemSimInitKey(){} ~ChemSimInitKey() {} };

///Key for adding and removing reactions
class ChemSimReactionKey {friend class CCylinder;
                          friend class Compartment;
                          friend class CompartmentGrid;
#ifdef TESTING
                          public:
#endif //TESTING
                          ChemSimReactionKey(){}; ~ChemSimReactionKey(){}; };
    
///Key for running the ChemSim
class ChemSimRunKey {friend class CController;
#ifdef TESTING
                     public:
#endif //TESTING
                     ChemSimRunKey(){} ~ChemSimRunKey() {} };
                    

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
    static void setInstance(ChemSimInitKey k, ChemSimImpl *csi);
    
    /// Copying is not allowed
    ChemSim(const ChemSim &rhs) = delete;

    /// Assignment is not allowed
    ChemSim& operator=(ChemSim &rhs) = delete;
    
    /// After all initial reactions have been added via addReaction(...) method, invoke initialize() prior to invoking run() 
    static void initialize(ChemSimInitKey k);
    
    /// Add Reaction *r to the chemical network which needs to be simulated
    static void addReaction(ChemSimReactionKey k, ReactionBase *r);
    
    /// Remove Reaction *r from the simulated chemical network 
    static void removeReaction(ChemSimReactionKey k, ReactionBase *r);
    
    /// Run the chemical dynamics for a specific number of steps
    static bool run(ChemSimRunKey k, int steps);
    
    /// Mainly used for debugging: print chemical reactions in the network at this moment
    static void printReactions();
private:

    static ChemSimImpl* _pimpl; ///< Store a pointer to a specific implementation of stochastic chemical kinetics; no ownership
    ChemSim() {};
};
#endif
