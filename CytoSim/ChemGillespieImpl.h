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

namespace chem {
    
    class RNodeGillespie;
    class ChemGillespieImpl;
    
    /// RNodeGillespie stands for Reaction Node Next Reaction Method.
    
    /*! RNodeGillespie manages a single chemical reaction within the Gillespie algorithm. It has a pointer to the PQ element
     *  containing the Reaction via a handle_t object (and hence can modify both the corresponding PQNode, such as PQNode's tau
     *  or the underlying Reaction instance). RNodeGillespie can recompute tau if needed and has auxiliary methods for computing
     *  reaction's propensity. When the propensity drops to zero, the RNodeGillespie can execute the passivateReaction() method.
     *  Alternatively, passivated RNodeGillespie can be activated via activateReaction(). The main part of the Gillespie algoritm is
     *  implemented in the makeStep() method.
     */
    class RNodeGillespie : public RNode {
    public:
        /// Ctor:
        /// @param *r is the Reaction object corresponding to this RNodeGillespie
        /// @param &chem_Gillespie is a refernce to ChemGillespieImpl object, which does the overall management of the Gillespie scheme (e.g. it
        /// gives acces to the heap itself, random distribution generators, etc.)
        RNodeGillespie(Reaction *r, ChemGillespieImpl &chem_Gillespie);
        
        /// Copying is not allowed
        RNodeGillespie(const RNodeGillespie& rhs) = delete;
        
        /// Assignment is not allowed
        RNodeGillespie& operator=(RNodeGillespie &rhs) = delete;
        
        /// Dtor: 1) Erases the corresponding PQNode element in the heap via the handle; 2) The RNode pointer of the
        /// tracked Reaction object is set to nullptr
        virtual ~RNodeGillespie();
                
        /// Returns a pointer to the Reaction which corresponds to this RNodeGillespie.
        Reaction* getReaction() const {return _react;};
        
        /// Return the currently stored propensity, "a", for this Reaction.
        /// @note The propensity is not recomputed in this method, so it potentially can be out of sync.
        double getPropensity() const {return _a;}
        
        void setA(double a) {_a=a;}
        
        /// Return the propensity, "a", associated with the penultimate step of this Reaction.
        /// @note The propensity is not recomputed in this method, so it potentially can be out of sync.
        double getPenultStepPropensity() const {return _a_prev;}

        void setPenultA(double a_prev) {_a_prev=a_prev;}
        
        /// (Re)Compute and return the propensity associated with this Reaction.
        double reComputePropensity() {
            _a_prev=_a;
            _a=_react->computePropensity();
            return _a;
        }
        
        void reset() {
            _a_prev = 0;
            _a = _react->computePropensity();
        }
        
        /// Compute and return the product of reactant copy numbers. This method can be used to quickly check
        /// whether this reaction needs to be passivated, if the returned result is zero.
        int getProductOfReactants () {return _react->getProductOfReactants ();};
        
        /// This method calls the corresponding Reaction::makeStep() method of the underyling Reaction object
        void makeStep() {_react->makeStep();}
        
        /// When this method is called, a new tau is computed and the corresponding PQNode element is updated in the heap.
        /// @note This method does not activate the Reaction itself, but instead only deals with the activation
        ///       process related to the corresponding PQNode element.
        virtual void activateReaction();
        
        /// When this method is called, reaction's tau is set to infinity, the propensity is set to 0, and
        /// the corresponding PQNode element is updated in the heap.
        /// @note This method does not passivate the Reaction itself, but instead only deals with the activation
        ///        process related to the corresponding PQNode element.
        virtual void passivateReaction();
        
        /// Return true if the Reaction is currently passivated
        bool isPassivated() const {return _react->isPassivated();}
        
        /// Print information about this RNodeGillespie such as tau, a and the Reaction which this RNodeGillespie tracks.
        void printSelf() const;
        
        /// Print the RNode objects which are dependents of this RNode (via the tracked Reaction object dependencies)
        void printDependents() const;
    private:
        ChemGillespieImpl &_chem_Gillespie; ///< A reference to the ChemGillespieImpl which containts the heap, random number generators, etc.
        Reaction *_react; ///< The pointer to the associated Reaction object. The corresponding memory is not managed by RNodeGillespie.
        double _a; ///< The propensity associated with the Reaction. It may be outdated and may need to be recomputed if needed.
        double _a_prev; ///< The propensity associated with the penultimate step of this Reaction.
    };
    
    
    
    /// ChemGillespieImpl stands for Chemical Next Reaction Method Implementation.
    
    /*! ChemGillespieImpl manages the Gillespie algorithm at the level of the network of reactions. In particular, this class contains
     *  the Gillespie heap and the exponential random number generator. Reaction objects can be added and removed from the
     *  ChemGillespieImpl instance.
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
        /// When the RNodeGillespie objects are destructed, they in turn destruct the corresponding PQNode element, setting the RNode
        /// pointer of the Reaction object to null. At the end, the heap itself will go out of scope.
        virtual ~ChemGillespieImpl();
        
        /// Return the number of reactions in the network.
        size_t getSize() const {return _n_reacts;}
        
        /// Return the current global time (as defined in the Gillespie algorithm)
        double getTime() const {return _t;}
        
        /// Sets global time to 0.0
        void resetTime() {_t=0.0; syncGlobalTime(); }
        
        /// Sets global time variable to ChemGillespieImpl's global time
        void syncGlobalTime() {global_time=_t;}
                
        /// Add Reaction *r to the network
        virtual void addReaction(Reaction *r);
        
        /// Remove Reaction *r from the network
        virtual void removeReaction(Reaction *r);
        
        double computeTotalA();
        
        /// A pure function (without sideeffects), which returns a random time tau, drawn from the exponential distribution,
        /// with the propensity given by a.
        double generateTau(double a);
        
        double generateUniform();
        
        /// This function iterates over all RNodeGillespie objects in the network, generating new tau-s for each case and
        /// subsequently updating the heap. It needs to be called before calling run(...).
        /// @note If somewhere in the middle of simulaiton initialize() is called, it will be analogous to starting
        /// the simulation from scratch, except with the Species copy numbers given at that moment in time. The global time
        /// is reset to zero again.
        virtual void initialize();
        
        /// This method runs the Gillespie algorithm for the given number of steps.
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
        
        void activateReaction(Reaction *r);
        
        void passivateReaction(Reaction *r);
        
        /// Prints all RNodes in the reaction network
        virtual void printReactions() const;
        
    private:
        /// This is a somewhat complex subroutine which implements the main part of the Gibson-Bruck Gillespie algoritm. See the
        /// implementation for details. After this method returns, roughly the following will have happened:
        /// 1) The Reaction corresponding to the lowest tau RNodeGillespie is executed and the corresponding Species copy numbers are changed
        /// 2) A new tau is computed from this Reaction and the corresponding PQNode element is updated in the heap
        /// 3) The other affected Reaction objects are found, their taus are recomputed and corresponding PQNode elements are
        ///    updated in the heap.
        /// 4) For the Reaction and associated Species signals are emitted, if these objects broadcast signals upon change.
        /// Returns true if successful, and false if the heap is exchausted and there no more reactions to fire
        bool makeStep();
    private:
        std::unordered_map<Reaction*, std::unique_ptr<RNodeGillespie>> _map_rnodes; ///< The database of RNodeGillespie objects, representing the reaction network
            std::unordered_map<Reaction*, std::unique_ptr<RNodeGillespie>> _map_rnodes_inactive; ///< The database of RNodeGillespie objects, representing the passivated part of the reaction network
        std::mt19937 _eng; ///< Random number generator
        std::exponential_distribution<double> _exp_distr; ///< Adaptor for the exponential distribution
        std::uniform_real_distribution<double> _uniform_distr;
        double _t; ///< global time
        double _a_total; 
        size_t _n_reacts; ///< number of reactions in the network
    };
    
} // end of namespace

#endif
