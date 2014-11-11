//
//  ChemSimImpl.h
//  CytoSim
//
//  Created by Garegin Papoian on 7/14/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_ChemSimImpl_h
#define CytoSim_ChemSimImpl_h

#include "common.h"

///FORWARD DECLARATIONS
class ReactionBase;

/// ChemSimImpl is an abstract base class for algorithms that run stochastic chemical kinetics.  

/*! Specific stochastic kinetics algorithm classes should inherit from ChemSimImpl. A user will 
 *  then attach the corresponding algorithm to ChemSim via the algoritm's base class ChemSimImpl.
 */
class ChemSimImpl {
public:
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as potentially throwing,
    /// which in turn disables move operations by the STL containers. This behaviour is a gcc bug
    /// (as of gcc 4.703), and will presumbaly be fixed in the future.
    virtual ~ChemSimImpl() noexcept {};
    
    /// After all initial reactions have been added via addReaction(...) method, invoke initialize() prior to invoking run() 
    virtual void initialize() = 0;
    
    /// Add Reaction *r to the chemical network which needs to be simulated
    virtual void addReaction(ReactionBase *r) = 0;
    
    /// Remove Reaction *r from the simulated chemical network 
    virtual void removeReaction(ReactionBase *r) = 0;
    
    /// Run the chemical dynamics for a specific number of steps
    virtual bool run(int steps) = 0;
    
    /// Mainly used for debugging: print chemical reactions in the network at this moment
    virtual void printReactions() const = 0;
};

#endif
