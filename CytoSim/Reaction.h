//
//  Reaction.h
//  CytoSim-Experimenting
//
//  Created by Garegin Papoian on 4/19/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_Experimenting_Reaction_h
#define CytoSim_Experimenting_Reaction_h

#include <iostream>
#include <array>
#include <algorithm>
#include <utility>
#include <unordered_set>
#include <cmath>
#include <initializer_list>

#include <boost/signals2/signal.hpp>
#include <boost/signals2/connection.hpp>
#include <boost/signals2/shared_connection_block.hpp>

#include "common.h"
#include "Species.h"

namespace chem {

    class RNode;
    class Composite;
    class ReactionBase;
    class SpeciesPtrContainerVector;

    
/// This is a ReactionBase signal object that will be called by ChemSignal, usually when requested by the 
/// ReactionBase simulation algorithm
typedef boost::signals2::signal<void (ReactionBase *)> ReactionEventSignal;


/// ReactionBase class represents simple chemical ReactionBases of the form A + B -> C.  

/*! ReactionBase is defined in terms of reactant and product RSpecies. In the current implementation, 
 *  the stoichiometric coefficents can only be one: i.e. A + B -> C + D is allowed, but A + 2B -> C + D is not allowed. 
 *  Almost all chemical ReactionBases are either unimolecular or bimolecular, so this restriction should not be too
 *  burdensom. Also, the ReactionBase indicates a forward process only. For processes in both directions, e.g. A <-> B, 
 *  two ReactionBases objects need to be defined, corresponding to A->B and B->A. 
 *
 *  A ReactionBase tracks other ReactionBase objects that are affected if this ReactionBase is executed. A ReactionBase may be set 
 *  up such that it "signals" when a ReactionBase event happens, in which case the corresponding callbacks are called.
 *
 */
        
class ReactionBase {
protected:
    std::vector<ReactionBase*> _dependents; ///< Pointers to ReactionBase objects that depend on this ReactionBase being executed
    RNode* _rnode; ///< A pointer to an RNode object which is used to implement a Gillespie-like algorithm (e.g. NRM)
    float _rate; ///< the rate for this ReactionBase
    // Composite *_parent;
#ifdef REACTION_SIGNALING
    ReactionEventSignal* _signal; ///< Can be used to broadcast a signal associated with this ReactionBase (usuall when a single step of this ReactionBase occurs)
#endif
#if defined TRACK_ZERO_COPY_N || defined TRACK_UPPER_COPY_N
    bool _passivated; ///< Indicates whether the ReactionBase is currently passivated
#endif 
    
public:
//    /// Default Constructor produces a ReactionBase with no reactants or products, zero rate, etc.
//    ReactionBase () :  _rnode(nullptr), _rate(0.0), _m(0), _is_signaling (false) {}  
//    

    /// The main constructor:
    /// @param Species that are reactants and products are put together into a single list (starting from reactants)
    /// @param M - number of reactants
    /// @param N - number of products
    /// @param rate - the rate constant for this ReactionBase
    ReactionBase (float rate);
        
    ReactionBase (const ReactionBase &rb) = delete; // no copying (including all derived classes)
    ReactionBase& operator=(ReactionBase &rb) = delete;  // no assignment (including all derived classes)

    /// Destructor
    virtual ~ReactionBase() noexcept
    {
#ifdef REACTION_SIGNALING
        if(_signal!=nullptr)
            delete _signal;
#endif
    }

    ReactionBase* clone(const SpeciesPtrContainerVector &spcv) {
        return cloneImpl(spcv);
    }
    
    virtual ReactionBase* cloneImpl(const SpeciesPtrContainerVector &spcv) = 0;
    
    /// Sets the ReactionBase rate to the parameter "rate" 
    void setRate(float rate) {_rate=rate;}
    
    /// Sets the RNode pointer associated with this ReactionBase to rhs. Usually is called only by the 
    /// Gillespie-like algorithms.
    void setRnode(RNode *rhs) {_rnode=rhs;} 

    /// Returns the rate associated with this ReactionBase
    float getRate() const {return _rate;}
    
    /// Returns a pointer to the RNode associated with this ReactionBase
    RNode* getRnode() const {return _rnode;} 
    
    /// Returns the number of reactant RSpecies
    unsigned short getM() const {return getMImpl();}
    
    virtual unsigned short getMImpl() const = 0;
    
    /// Returns the number of product RSpecies
    unsigned short getN() const {return getNImpl();}
    
    virtual unsigned short getNImpl() const = 0;
        
//    Composite* getParent() {return _parent;}
//    
//    void setParent (Composite *other) {_parent=other;}
//    
//    bool hasParent() const {return _parent!=nullptr? true : false;}
//    
//    Composite* getRoot();
    
    Composite* getParent() {return nullptr;}
    
    void setParent (Composite *other) {}
    
    bool hasParent() const {return false;}
    
    Composite* getRoot() {return nullptr;}
    
    /// Computes the product of the copy number of all reactant RSpecies.
    /// Can be used to quickly determine whether this ReactionBase should be allowed to activate - if one of the
    /// reactants has a copy number equal to zero, then zero is returned, indicating that
    /// this ReactionBase should not be (yet) activated.
    int getProductOfReactants () const {return getProductOfReactantsImpl();}
    
    virtual int getProductOfReactantsImpl() const = 0;
    
    /// Computes the product of the copy number of all product RSpecies minus maxium allowed copy number.
    /// Can be used to quickly determine whether this ReactionBase should be allowed to activate - if one of the
    /// products has a copy number equal to the maximum allowed, then zero is returned, indicating that
    /// this ReactionBase should not be (yet) activated.
    int getProductOfProducts () const {return getProductOfProductsImpl();}
    
    virtual int getProductOfProductsImpl() const = 0;
    
    /// Return true if the ReactionBase is currently passivated
#if defined TRACK_ZERO_COPY_N || defined  TRACK_UPPER_COPY_N
    bool isPassivated() const {return _passivated;}
#else
    bool isPassivated() const {return false;}
#endif
    
    bool containsSpecies(Species *s) const {return containsSpeciesImpl(s);}
    
    virtual bool containsSpeciesImpl(Species *s) const = 0;
    
#ifdef REACTION_SIGNALING
    /// Return true if this RSpecies emits signals on copy number change
    bool isSignaling () const {return _signal!=nullptr;}
    
    /// Set the signaling behavior of this ReactionBase
    void startSignaling ();
    
    /// Destroy the signal associated with this ReactionBase; all associated slots will be destroyed
    /// @note To start signaling again, startSignaling() needs to be called
    void stopSignaling ();

    /// Connect the callback, react_callback to a signal corresponding to ReactionBase *r.
    /// @param std::function<void (ReactionBase *)> const &react_callback - a function object to be called (a slot)
    /// @param int priority - lower priority slots will be called first. Default is 5 Do not use priorities 1 and 2 
    ///                       unless absolutely essential.
    /// @return a connection object which can be used to later disconnect this particular slot or temporarily block it
    boost::signals2::connection connect(std::function<void (ReactionBase *)> const &react_callback, int priority=5);
    
    /// Broadcasts signal indicating that the ReactionBase event has taken place 
    /// This method is only called by the code which runs the chemical dynamics (i.e. Gillespie-like algorithm)
    void emitSignal() {
        if(isSignaling())
            (*_signal)(this);
    }
#endif
    
    /// Return a const reference to the vector of dependent ReactionBases
    /// @note One can obtain two different lists of affected ReactionBases:
    /// 1) via getAffectedReactionBases(), where the copy numbers do influence the 
    /// dependencies, and 2) via dependents(), where dependencies stop being counted 
    /// if the copy numbers of reactant species drop to 0. 
    const std::vector<ReactionBase*>& dependents() {return _dependents;}
            
    /// Returns true if two ReactionBase objects are equal.
    /// Two ReactionBase objects are equal if each of their reactants and products are equal
    friend bool operator==(const ReactionBase& a, const ReactionBase& b)
    {
        if(typeid(a) != typeid(b))
            return false;
        return a.is_equal(b); // Invoke virtual is_equal via derived subclass of a (which should be the same as b)
    }
    
    virtual bool is_equal(const ReactionBase& b) const = 0;
        
    /// Return true if two ReactionBase are not equal.
    /// @see operator ==(const ReactionBase& a, const ReactionBase& b) above
    friend bool operator !=(const ReactionBase& a, const ReactionBase& b){
        return !(a==b);
    }
    
    /// Fire the ReactionBase - make a single step, where reactant RSpecies copy numbers are
    /// decreased by one, and the product RSpecies copy numbers are increased by one.
    /// @note This method does not send a ReactionBase event Signal. The latter is usually done
    /// from within a Gillespie-like algorithm. 
    void makeStep() {makeStepImpl();}
    
    virtual void makeStepImpl() = 0;
    
    /// Compute the ReactionBase propensity that is needed by a Gillespie like algorithm:
    /// rate*reactant_1.getN()*reactant_2.getN()...
    float computePropensity() const {return computePropensityImpl();}

    virtual float computePropensityImpl() const = 0;
    
    /// Usually is applied to ReactionBase objects with propensity of 0 (e.g. when one of the 
    /// copy numbers of reactants has dropped to 0. This method call notifies all other 
    /// ReactionBase objects that may affect this ReactionBase to stop tracking this ReactionBase. Eventually, 
    /// activateReaction() may be called to restart tracking, if the propensity stops being 0.
    void passivateReaction() {passivateReactionImpl();}
    
    virtual void passivateReactionImpl() = 0;
    
    /// Requests that ReactionBase objects that may affect this Reaction to start tracking it, which can be
    /// used to follow ReactionBase objects whose propensities change upong firing of some ReactionBase. This
    /// request is acted upond unconditionally.
    void activateReactionUnconditional() {activateReactionUnconditionalImpl();}
    
    virtual void activateReactionUnconditionalImpl() = 0;
    
    /// Requests that Reaction objects that may affect this Reaction to start tracking it, which can be
    /// used to follow Reaction objects whose propensities change upong firing of some Reaction. This
    /// request will be ignored if the Reaction's propensity is still zero.
    void activateReaction() {
#ifdef TRACK_ZERO_COPY_N
        if(getProductOfReactants()==0) // One of the reactants is still at zero copy n, no need to activate yet...
            return;
#endif
#ifdef TRACK_UPPER_COPY_N
        if(getProductOfProducts()==0) // One of the products is at the maximum allowed copy number, no need to activate yet...
            return;
#endif
        activateReactionUnconditional();
    }
    
    /// Print the ReactionBases that are affacted by this ReactionBase being fired
    void printDependents() ;
    
    /// Return the list of ReactionBase objects that are affected when this ReactionBase is fired
    /// @return a vector of pointers to the affected ReactionBase objects
    /// @note This method is "expensive" because it computes from scratch the dependencies. Importantly, 
    /// the copy numbers of molecules do not influence the result of this function. \sa dependents()
    virtual std::vector<ReactionBase*> getAffectedReactions() = 0;
    
    /// Request that the ReactionBase *r adds this ReactionBase to its list of dependents which it affects. 
    void registerNewDependent(ReactionBase *r);
    
    /// Request that the ReactionBase *r removes this ReactionBase from its list of dependents which it affects. 
    /// This is usually requested when the ReactionBase propensity drops to zero (i.e. via passivateReactionBase()).
    void unregisterDependent(ReactionBase *r);

    virtual void printToStream(std::ostream& os) const = 0;
    
    /// Print self into an iostream
    friend std::ostream& operator<<(std::ostream& os, const ReactionBase& rr)
    {
        rr.printToStream(os);
        return os;
    }
};
    
template <unsigned short M, unsigned short N>
    class Reaction : public ReactionBase {
    private:
        std::array<RSpecies*, M+N> _rspecies;
    public:
        Reaction(const std::vector<Species*> &species, float rate = 0.0) : ReactionBase(rate)
        {
            initializeSpecies(species);
        }
        
        Reaction (const Reaction &rb) = delete; // no copying (including all derived classes)
        Reaction& operator=(Reaction &rb) = delete;  // no assignment (including all derived classes)

#ifdef BOOST_MEM_POOL
        /// Advanced memory management
        void* operator new(std::size_t size);
        
        void operator delete(void* ptr) noexcept;
#endif
        void initializeSpecies(const std::vector<Species*> &species);

        /// Destructor
        virtual ~Reaction() noexcept
        {
        for(auto i=0U; i<M; ++i)
            _rspecies[i]->removeAsReactant(this);
        for(auto i=M; i<(M+N); ++i)
            _rspecies[i]->removeAsProduct(this);

        }
        std::array<RSpecies*, M+N>& rspecies() {return _rspecies;}
        const std::array<RSpecies*, M+N>& rspecies() const {return _rspecies;}
        
        virtual std::vector<ReactionBase*> getAffectedReactions() override
        {
            std::unordered_set<ReactionBase*> rxns;
            for(auto s : _rspecies){
                //                std::cout << "std::vector<ReactionBase*> getAffectedReactions(): " << *s << std::endl;
                rxns.insert(s->beginReactantReactions(),s->endReactantReactions());
            }
            //        std::sort(rxns.begin(),rxns.end());
            //        rxns.erase(std::unique(rxns.begin(),rxns.end()), rxns.end());
            rxns.erase(this);
            return std::vector<ReactionBase*>(rxns.begin(),rxns.end());
        }
        
    protected:
        virtual unsigned short getMImpl() const override {return M;};
        virtual unsigned short getNImpl() const override {return N;};
        
        virtual bool is_equal(const ReactionBase& other) const override
        {
            const Reaction<M,N> *a = this;
            const Reaction<M,N> *b = static_cast<const Reaction<M,N>*>(&other);
            auto it_pair = std::mismatch(a->_rspecies.begin(),a->_rspecies.end(),b->_rspecies.begin(),
                                         [](RSpecies* A, RSpecies* B){return A->getSpecies()==B->getSpecies();});
            if(it_pair.first==a->_rspecies.end())
                return true;
            return false;
        }
        
        virtual int getProductOfReactantsImpl() const override
        {
            int prod = 1;
            for(auto i=0U; i<M; ++i)
                prod*=_rspecies[i]->getN();
            return prod;
            
        }

        virtual float computePropensityImpl() const override
        {
#ifdef TRACK_UPPER_COPY_N
            if(this->Reaction<M,N>::getProductOfProductsImpl()==0){
                //            std::cout << "Reaction::computePropensity() for the Reaction, " << (*this)
                //            << " will return 0.0";
                return float(0.0);
            }
#endif
//            return _rate*Reaction<M,N>::getProductOfReactantsImpl();
            int prod = 1;
            for(auto i=0U; i<M; ++i)
                prod*=_rspecies[i]->getN();
            return _rate*prod;
        }
        
        virtual int getProductOfProductsImpl() const override
        {
#ifdef TRACK_UPPER_COPY_N
            int prod = 1;
            for(auto i=M; i<(M+N); ++i){
                prod*=_rspecies[i]->getN()-_rspecies[i]->getUpperLimitForN();
                //            std::cout << "getProductOfProducts(): " << (*this) << (*sit)->getN() << " " << (*sit)->getUpperLimitForN() << " " << ((*sit)->getN()-(*sit)->getUpperLimitForN()) << std::endl;
            }
            return prod;
#else
            return 1;
#endif
        }
    
        virtual bool containsSpeciesImpl(Species *s) const override
        {
            auto it = std::find_if(_rspecies.begin(), _rspecies.end(),
                                   [s](const RSpecies *rs){return (&rs->getSpecies())==s;});
            return it!=_rspecies.end();
            
        }
    
        virtual void makeStepImpl() override
        {
            for(auto i=0U; i<M; ++i)
                _rspecies[i]->down();
            for(auto i=M; i<(M+N); ++i)
                _rspecies[i]->up();
        }

        virtual void activateReactionUnconditionalImpl() override;
        virtual void passivateReactionImpl() override;
        
        virtual void printToStream(std::ostream& os) const override
        {
            unsigned char i=0;
            auto sit = _rspecies.cbegin();
            auto send = sit+M;
            for (; sit!=send; ++sit)
            {
                os << (*sit)->getFullName() << "{" << (*sit)->getN() << "}";
                if(i<M-1)
                os << " + ";
                ++i;
            }
            os << " ---> ";
            i=0;
            for (auto sit = send; sit!=_rspecies.cend(); ++sit)
            {
                os << (*sit)->getFullName() << "{" << (*sit)->getN() << "}";
                if(i<((N-M)-1))
                    os << " + ";
                ++i;
            }
            os << ", " << "curr_rate = " << getRate() << ", a=" << computePropensity() << ", ReactionBase ptr=" << this << "\n";
        }

        virtual Reaction<M,N>* cloneImpl(const SpeciesPtrContainerVector &spcv) override;
    };
                

} // end of namespace
#endif
