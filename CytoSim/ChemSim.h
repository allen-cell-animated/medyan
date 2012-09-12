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
#include "Reaction.h"

namespace chem {
    
class ChemSimImpl;

/// ChemSim is used to manage running a network of chemical reactions.  

/*! ChemSim implements a Strategy pattern, allowing custom algorithms for running stochastic chemical simulations. 
 *  After the specific algorithm is chosen and ChemSim is instantiated, ChemSim can be used to manage simulations, through 
 *  such methods as run(steps) etc. Here is an example:
 *  @code
        SpeciesBulk A1("A1",  25);
        SpeciesBulk A2("A2", 25);
        Reaction r1 = { {&A1,&A2}, 1, 1, 100.0 };
        ChemNRMImpl chem_nrm_impl;
        ChemSim chem(&chem_nrm_impl);
        chem.addReaction(&r1);
        chem.initialize();
        chem.printReactions();
        chem.run(100);   
        chem.printReactions();
 * @endcode
 */
class ChemSim {
public:
    /// Constructor 
    /// @param ChemSimImpl *csi is a pointer the concrete implementation of stochastic simulation algorithm.
    /// @note ChemSim simply stores the csi pointer but does not manage its memory. Make sure that csi is always 
    /// a valid pointer while ChemSim is used.
    ChemSim(ChemSimImpl *csi);

    /// Copying is not allowed
    ChemSim(const ChemSim &rhs) = delete;

    /// Assignment is not allowed
    ChemSim& operator=(ChemSim &rhs) = delete;
    
    /// After all initial reactions have been added via addReaction(...) method, invoke initialize() prior to invoking run() 
    void initialize();
    
    /// Add Reaction *r to the chemical network which needs to be simulated
    void addReaction(ReactionBase *r);
    
    /// Remove Reaction *r from the simulated chemical network 
    void removeReaction(ReactionBase *r);
    
    /// Run the chemical dynamics for a specific number of steps
    bool run(int steps);
    
    /// Mainly used for debugging: print chemical reactions in the network at this moment
    void printReactions() const;
private:
    ChemSimImpl* _pimpl; ///< Store a pointer to a specific implementation of stochastic chemical kinetics; no ownership
};

} // end of namespace
#endif
