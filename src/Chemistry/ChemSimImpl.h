
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

#ifndef MEDYAN_ChemSimImpl_h
#define MEDYAN_ChemSimImpl_h

#include "common.h"
#include "Histogram.h"
#include "DissipationTracker.h"
#include <fstream>

//FORWARD DECLARATIONS
class ReactionBase;

/// An abstract base class for algorithms that run stochastic chemical kinetics.  

/*! Specific stochastic kinetics algorithm classes should inherit from ChemSimImpl. 
 *  A user will then attach the corresponding algorithm to ChemSim via the algoritm's 
 *  base class ChemSimImpl.
 */
class ChemSimImpl {
public:
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~ChemSimImpl() noexcept {};
    
    /// After all initial reactions have been added via addReaction(...) method, invoke
    /// initialize() prior to invoking run()
    virtual void initialize() = 0;
    
    /// Add Reaction *r to the chemical network which needs to be simulated
    virtual void addReaction(ReactionBase *r) = 0;
    
    /// Remove Reaction *r from the simulated chemical network 
    virtual void removeReaction(ReactionBase *r) = 0;
    
    /// Run the chemical dynamics for a set amount of time
    virtual bool run(double time) = 0;
    
    /// Run the chemical dynamics for a set amount of reaction steps
    virtual bool runSteps(int steps) = 0;
    
    /// Mainly used for debugging: print chemical reactions in the network at
    /// this moment
    virtual void printReactions() const = 0;
    
    DissipationTracker * _dt = nullptr;
    
};

#endif
