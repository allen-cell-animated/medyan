//
//  ChemSimpleGillespieImpl.h
//  CytoSim
//
//  Created by Garegin Papoian on 8/16/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_ChemSimpleGillespieImpl_h
#define CytoSim_ChemSimpleGillespieImpl_h

#include <vector>
#include <random>

#include "utility.h"
#include "Reaction.h"
#include "ChemSimImpl.h"
#include "ChemRNode.h"

namespace chem {
    
    /// ChemSimpleGillespieImpl implements an unoptimized version of the Gillespie algorithm.
    
    /*! ChemSimpleGillespieImpl manages the Gillespie algorithm at the level of the network of reactions. In particular, this class contains
     *  the Gillespie heap and the exponential random number generator. Reaction objects can be added and removed from the
     *  ChemSimpleGillespieImpl instance.
     */
    class ChemSimpleGillespieImpl : public ChemSimImpl {
    public:
        /// Ctor: Seeds the random number generator, sets global time to 0.0 and the number of reactions to 0
        ChemSimpleGillespieImpl() :
        ChemSimImpl(), _eng(static_cast<unsigned long>(time(nullptr))), _exp_distr(0.0), _uniform_distr() {
            resetTime();
        }
        
        /// Copying is not allowed
        ChemSimpleGillespieImpl(const ChemSimpleGillespieImpl &rhs) = delete;
        
        /// Assignment is not allowed
        ChemSimpleGillespieImpl& operator=(ChemSimpleGillespieImpl &rhs) = delete;
        
        ///Dtor: The reaction network is cleared. The RNodeGillespie objects will be destructed, but Reaction objects will stay intact.
        /// When the RNodeGillespie objects are destructed, they in turn destruct the corresponding PQNode element, setting the RNode
        /// pointer of the Reaction object to null. At the end, the heap itself will go out of scope.
        virtual ~ChemSimpleGillespieImpl();
        
        /// Return the number of reactions in the network.
        size_t getSize() const {return _reactions.size();}
        
        /// Return the current global time (as defined in the Gillespie algorithm)
        double getTime() const {return _t;}
        
        /// Sets global time to 0.0
        void resetTime() {_t=0.0; syncGlobalTime(); }
        
        /// Sets global time variable to ChemSimpleGillespieImpl's global time
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
        void initialize();
        
        /// This method runs the Gillespie algorithm for the given number of steps.
        virtual bool run(int steps) {
            for(int i=0; i<steps; ++i){
                bool success = makeStep();
                if(!success)
                    return false;
                //            if(i%1000000==0)
                //                std::cout << "ChemSimpleGillespieImpl::run(): i=" << i << std::endl;
            }
            return true;
        }
        
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
        std::vector<Reaction*> _reactions; ///< The database of Reaction objects, representing the reaction network
        std::mt19937 _eng; ///< Random number generator
        std::exponential_distribution<double> _exp_distr; ///< Adaptor for the exponential distribution
        std::uniform_real_distribution<double> _uniform_distr;
        double _t; ///< global time
    };
    
} // end of namespace

#endif
