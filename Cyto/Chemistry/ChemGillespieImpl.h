//
//  ChemGillespieImpl.h
//  CytoSim
//
//  Created by Garegin Papoian on 8/11/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_ChemGillespieImpl_h
#define CytoSim_ChemGillespieImpl_h

#include <vector>
#include <random>

#include "utility.h"
#include "Reaction.h"
#include "ChemSimImpl.h"
#include "ChemRNode.h"

class RNodeGillespie;
class ChemGillespieImpl;
    
/// RNodeGillespie is used by ChemGillespieImpl to implement the cached version of the Gillespie algorithm.

/*! RNodeGillespie manages a single chemical reaction within the Gillespie algorithm. When the propensity drops to zero, the RNodeGillespie can execute the passivateReaction() method. Alternatively, passivated RNodeGillespie can be activated via activateReaction(). The main part of the Gillespie algoritm is implemented in the makeStep() method.
 */
class RNodeGillespie : public RNode {
public:
    /// Ctor:
    /// @param *r is the Reaction object corresponding to this RNodeGillespie
    /// @param &chem_Gillespie is a refernce to ChemGillespieImpl object, which does the overall management of the Gillespie scheme (e.g.   random distribution generators, etc.)
    RNodeGillespie(ReactionBase *r, ChemGillespieImpl &chem_Gillespie);
    
    /// Copying is not allowed
    RNodeGillespie(const RNodeGillespie& rhs) = delete;
    
    /// Assignment is not allowed
    RNodeGillespie& operator=(RNodeGillespie &rhs) = delete;
    
    /// Dtor: The RNode pointer of the tracked Reaction object is set to nullptr
    virtual ~RNodeGillespie();
            
    /// Returns a pointer to the Reaction which corresponds to this RNodeGillespie.
    ReactionBase* getReaction() const {return _react;};
    
    /// Return the currently stored propensity, "a", for this Reaction.
    /// @note The propensity is not recomputed in this method, so it potentially can be out of sync.
    double getPropensity() const {return _a;}
    
    /// Set the propensity, "a", for this Reaction.
    void setA(double a) {_a=a;}
    
    /// Return the propensity, "a", associated with the penultimate step of this Reaction.
    /// @note The propensity is not recomputed in this method, so it potentially can be out of sync.
    double getPenultStepPropensity() const {return _a_prev;}

    /// Set the propensity, "a", associated with the penultimate step of this Reaction.
    void setPenultA(double a_prev) {_a_prev=a_prev;}
    
    /// (Re)Compute and return the propensity associated with this Reaction.
    /// Remembers the penultimate step propensity as well
    double reComputePropensity() {
        _a_prev=_a;
        _a=_react->computePropensity();
        return _a;
    }
    
    /// Set the the penultimate step propensity to zero and compute
    /// the current propensity.
    void reset() {
        _a_prev = 0;
        _a = _react->computePropensity();
    }
    
    /// This method calls the corresponding Reaction::makeStep() method of the underyling Reaction object
    void makeStep() {_react->makeStep();}
    
    /// Forwards the call to the similarly named method of ChemGillespieImpl
    virtual void activateReaction();
    
    /// Forwards the call to the similarly named method of ChemGillespieImpl
    virtual void passivateReaction();
    
    /// Forwards the call to the similarly named method of ChemGillespieImpl
    bool isPassivated() const {return _react->isPassivated();}
    
    /// Print information about this RNodeGillespie such as "a", "a_penult" and the Reaction which this RNodeGillespie tracks.
    void printSelf() const;
    
    /// Print the RNode objects which are dependents of this RNode (via the tracked Reaction object dependencies)
    void printDependents() const;
private:
    ChemGillespieImpl &_chem_Gillespie; ///< A reference to the ChemGillespieImpl which containts the heap, random number generators, etc.
    ReactionBase *_react; ///< The pointer to the associated Reaction object. The corresponding memory is not managed by RNodeGillespie.
    double _a; ///< The propensity associated with the Reaction. It may be outdated and may need to be recomputed if needed.
    double _a_prev; ///< The propensity associated with the penultimate step of this Reaction.
};



/// ChemGillespieImpl implements a slightly optimized version of the Gillespie algorithm.

/*! ChemGillespieImpl manages the Gillespie algorithm at the level of the network of reactions. Reaction objects can be added and removed from the ChemGillespieImpl instance. The propensities of all Reactions are cached, and they are recomputed only when the copy number of associated reactant species gets changed (due to the currently occurring Reaction). 
    @note The algorithm used by this class relies on tracking dependent 
    Reactions, so TRACK_DEPENDENTS must be defined. The algorithm can work 
    both when TRACK_ZERO_COPY_N and TRACK_UPPER_COPY_N are either defined or 
    undefined. In the former case, Reaction objects with zero propensities 
    stop being treated as dependents until their propensities become nonzero again.
 */
class ChemGillespieImpl : public ChemSimImpl {
public:
    /// Ctor: Seeds the random number generator, sets global time to 0.0 and the number of reactions to 0
    ChemGillespieImpl() :
    ChemSimImpl(), _eng(static_cast<unsigned long>(time(nullptr))), _exp_distr(0.0), _uniform_distr(), _a_total(0),_n_reacts(0) {
        resetTime();
    }
    
    /// Copying is not allowed
    ChemGillespieImpl(const ChemGillespieImpl &rhs) = delete;
    
    /// Assignment is not allowed
    ChemGillespieImpl& operator=(ChemGillespieImpl &rhs) = delete;
    
    ///Dtor: The reaction network is cleared. The RNodeGillespie objects will be destructed, but Reaction objects will stay intact.
    virtual ~ChemGillespieImpl();
    
    /// Return the number of reactions in the network.
    size_t getSize() const {return _map_rnodes.size();}
    
    /// Return the current global time (as defined in the Gillespie algorithm)
    double getTime() const {return _t;}
    
    /// Sets the global time to 0.0
    void resetTime() {_t=0.0; syncGlobalTime(); }
    
    /// Sets global time variable to ChemGillespieImpl's global time
    void syncGlobalTime() {global_time=_t;}
            
    /// Add ReactionBase *r to the network
    virtual void addReaction(ReactionBase *r);
    
    /// Remove ReactionBase *r from the network
    virtual void removeReaction(ReactionBase *r);
    
    /// Unconditionally compute the total propensity associated with the network.
    double computeTotalA();
    
    /// Returns a random time tau, drawn from the exponential distribution,
    /// with the propensity given by a.
    double generateTau(double a);
    
    /// Returns a random number between 0 and 1, drawn from the uniform distribution
   double generateUniform();
    
    /// This function iterates over all RNodeGillespie objects in the network, activating all Reaction objects and calling reset(). The total propentsity for the network is computed.
    /// @note This method needs to be called before calling run(...).
    /// @note If somewhere in the middle of simulaiton initialize() is called, it will be analogous to starting the simulation from scratch, except with the Species copy numbers given at that moment in time. The global time is reset to zero again.
    virtual void initialize();
    
    /// This method runs the Gillespie algorithm for the given number of steps.
    /// @return true if successful.
    virtual bool run(int steps) {
        for(int i=0; i<steps; ++i){
            bool success = makeStep();
            if(!success)
                return false;
            //            if(i%1000000==0)
            //                std::cout << "ChemGillespieImpl::run(): i=" << i << std::endl;
        }
        return true;
    }
    
    /// This method is used to track the change in the total propensity of the network as the previously passivated ReactionBase *r has become activated
    void activateReaction(ReactionBase *r);
    
    /// This method is used to track the change in the total propensity of the network as the ReactionBase *r has become passivated
    void passivateReaction(ReactionBase *r);
    
    /// Prints all RNodes in the reaction network
    virtual void printReactions() const;
    
private:

    /// This subroutine, along with with passivateReaction() and activateReaction() implements a cached version of the Gillespie algorith. Two random numbers get thrown, one for obtaining tau, the time to the next reaction event, and a uniform number to select which Reaction has occurred. Instead of computing the total propensity of the  network from scratch, the cached value is being modified as Reaction events occur.
    /// Returns true if successful, and false if the heap is exchausted and there no more reactions to fire
    bool makeStep();
private:
    std::unordered_map<ReactionBase*, std::unique_ptr<RNodeGillespie>> _map_rnodes; ///< The database of RNodeGillespie objects, representing the reaction network
    std::mt19937 _eng; ///< Random number generator
    std::exponential_distribution<double> _exp_distr; ///< Adaptor for the exponential distribution
    std::uniform_real_distribution<double> _uniform_distr;
    double _t; ///< global time
    double _a_total; 
    size_t _n_reacts; ///< number of reactions in the network
};



#endif
