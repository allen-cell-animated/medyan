
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_ChemSimpleGillespieImpl_h
#define MEDYAN_ChemSimpleGillespieImpl_h

#include <vector>
#include <random>

#include "common.h"

#include "Reaction.h"
#include "ChemSimImpl.h"
#include "ChemRNode.h"

/// Implements the simplest version of the Gillespie algorithm, without caching, etc.

/*! 
 *  ChemSimpleGillespieImpl manages the Gillespie algorithm at the level of the network
 *  of reactions. Reaction objects can be added and removed from the 
 *  ChemSimpleGillespieImpl instance.
 */
class ChemSimpleGillespieImpl : public ChemSimImpl {
public:
    /// Ctor: Seeds the random number generator, sets global time to 0.0
#ifdef DEBUGCONSTANTSEED
    ChemSimpleGillespieImpl() :
    ChemSimImpl(),
    _exp_distr(0.0), _uniform_distr() { resetTime(); }
#else
    ChemSimpleGillespieImpl() :
    ChemSimImpl(), _eng(rdtsc()),
    _exp_distr(0.0), _uniform_distr() { resetTime(); }
#endif
    
    /// Copying is not allowed
    ChemSimpleGillespieImpl(const ChemSimpleGillespieImpl &rhs) = delete;
    
    /// Assignment is not allowed
    ChemSimpleGillespieImpl& operator=(ChemSimpleGillespieImpl &rhs) = delete;
    
    ///Dtor: The reaction network is cleared.
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~ChemSimpleGillespieImpl() noexcept;
    
    /// Return the number of reactions in the network.
    size_t getSize() const {return _reactions.size();}
    
    /// Return the current global time (which should be the sum of all previously
    /// occurred tau-s of the Gillespie algorithm)
    double getTime() const {return _t;}
    
    /// Sets global time to 0.0
    void resetTime() {_t=0.0; syncGlobalTime(); }
    
    /// Sets global time variable to ChemSimpleGillespieImpl's global time
    void syncGlobalTime() {global_time=_t;}
    
    /// Add ReactionBase *r to the network
    virtual void addReaction(ReactionBase *r);
    
    /// Remove ReactionBase *r from the network
    virtual void removeReaction(ReactionBase *r);
    
    /// Compute the total propensity of the reaction network, by adding all individual
    /// reaction propensities
    double computeTotalA();
    
    /// A pure function (without sideeffects), which returns a random time tau, drawn
    /// from the exponential distribution, with the propensity given by a.
    double generateTau(double a);
    
    /// This function generates a random number between 0 and 1
    double generateUniform();
    
    /// This function needs to be called before calling run(...).
    /// @note If somewhere in the middle of simulaiton initialize() is called, it will
    /// be analogous to starting the simulation from scratch, except with the Species
    /// copy numbers given at that moment in time. The global time is reset to zero again.
    void initialize();
    
    /// This method runs the Gillespie algorithm for the given amount of time.
    /// @return true if successful.
    virtual bool run(double time) {
        
        double endTime = _t + time;
        
        while(_t < endTime) {
            bool success = makeStep();
            if(!success)
                return false;
        }
        return true;
    }
    
    /// This method runs the Gillespie algorithm for the given amount of reaction steps.
    /// @return true if successful.
    virtual bool runSteps(int steps) {
        
        for(int i = 0; i < steps; i++) {
            
            bool success = makeStep();
            if(!success)
                return false;
        }
        return true;
    }
    
    /// Prints all Reaction objects in the reaction network
    virtual void printReactions() const;
    
private:
    /// This subroutine implements the vanilla version of the Gillespie algorithm
    bool makeStep();
private:
    vector<ReactionBase*> _reactions; ///< The database of Reaction objects,
                                      ///< representing the reaction network
#ifndef DEBUGCONSTANTSEED
    mt19937 _eng; ///< Random number generator
#endif
    exponential_distribution<double> _exp_distr; ///< Adaptor for the exponential distribution
    uniform_real_distribution<double> _uniform_distr; ///< Adaptor for the uniform distribution
    double _t; ///< global time
};


#endif
