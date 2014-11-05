//
//  ReactionBase.h
//  CytoSim
//
//  Created by Garegin Papoian on 9/18/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef __CytoSim__ReactionBase__
#define __CytoSim__ReactionBase__

#include <iostream>

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

class RNode;
class Composite;
class ReactionBase;
class SpeciesPtrContainerVector;

///Enumeration for type of reaction
enum ReactionType {
    REGULAR, DIFFUSION, POLYMERIZATION, DEPOLYMERIZATION, BASICBINDING,
    LINKERBINDING, MOTORBINDING, BASICUNBINDING, LINKERUNBINDING, MOTORUNBINDING,
    MOTORWALKINGFORWARD, MOTORWALKINGBACKWARD,
    CREATION, DESTRUCTION
};

/// This is a ReactionBase signal object that may be called by a ReactionBase simulation algorithm
typedef boost::signals2::signal<void (ReactionBase *)> ReactionEventSignal;


/// ReactionBase class represents an abstract interface for simple chemical Reactions of the form A + B -> C.

/*! ReactionBase provides an interface for managing a chemical reaction. It is an abstract interface, so
 *  it cannot be directly instantiated, but other concrete classes may be derived from it. ReactionBase
 *  may have a composite object as a parent. A signaling interface may be used to make callbacks when
 *  some event, such as a single reaction step, has been executed.
 
 *  The ReactionBase indicates a forward process only. For processes in both directions, e.g. A <-> B,
 *  two ReactionBases objects need to be defined, corresponding to A->B and B->A.
 *
 *  A ReactionBase tracks other ReactionBase objects that are affected if this ReactionBase is executed. A ReactionBase may be set up such that it "signals" when a ReactionBase event happens, in which case the corresponding callbacks are called.
 *
 */

class ReactionBase {
protected:
    std::vector<ReactionBase*> _dependents; ///< Pointers to ReactionBase objects that depend on this ReactionBase being executed
    RNode* _rnode; ///< A pointer to an RNode object which is used to implement a Gillespie-like algorithm (e.g. NRM)
    Composite *_parent; ///< A pointer to a Composite object to which this Reaction belongs
    float _rate; ///< the rate for this ReactionBase
    float _rate_bare; ///< the bare rate for this ReactionBase (original rate)
#ifdef REACTION_SIGNALING
    std::shared_ptr<ReactionEventSignal> _signal; ///< Can be used to broadcast a signal associated with this ReactionBase (usuall when a single step of this ReactionBase occurs)
#endif
#if defined TRACK_ZERO_COPY_N || defined TRACK_UPPER_COPY_N
    bool _passivated; ///< Indicates whether the ReactionBase is currently passivated
#endif
    ReactionType _reactionType; ///< Reaction type enumeration
    bool _isProtoCompartment = false; ///Reaction is in proto compartment (Do not copy as a dependent, not in chemsim)
    
    ///For cross filament reactions only (should be moved to new class)
    float _rMin = 0.0;  ///< Reaction range
    float _rMax = 0.0;
    
    short _reactionID = -1; ///<Unique id of this reaction
    
public:
    /// The main constructor:
    /// @param rate - the rate constant for this ReactionBase
    ReactionBase (float rate, bool isProtoCompartment);
    
    /// No copying (including all derived classes)
    ReactionBase (const ReactionBase &rb) = delete;
    
    /// no assignment (including all derived classes)
    ReactionBase& operator=(ReactionBase &rb) = delete;
    
    /// Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as potentially throwing,
    /// which in turn disables move operations by the STL containers. This behaviour is a gcc bug
    /// (as of gcc 4.703), and will presumbaly be fixed in the future.
    virtual ~ReactionBase() noexcept {}
    
    /// Copy this reaction using SpeciesPtrContainerVector &spcv as a source of analogous Species.
    /// @return the cloned ReactionBase pointer.
    /// @note the receving object should take care of the memory management of the returned pointer
    ReactionBase* clone(const SpeciesPtrContainerVector &spcv) {
        return cloneImpl(spcv);
    }
    
    /// (Private) implementation of the clone() method to be elaborated in derived classes
    virtual ReactionBase* cloneImpl(const SpeciesPtrContainerVector &spcv) = 0;
    
    /// Returns a pointer to the first element of the std::array<RSpecies*>. This pointer can be used
    /// to iterate over RSpecies* if necessary (also by relying on getM() and size() to determine
    /// the iteration limits). The corresponding std::array<RSpecies*> is defined by the derived classes.
    virtual RSpecies** rspecies() = 0;
    
    ///Set reaction type
    void setReactionType(ReactionType rxnType) {_reactionType = rxnType;}
    
    ///Get reaction type
    ReactionType getReactionType() {return _reactionType;}
    
    /// Sets the ReactionBase rate to the parameter "rate"
    void setRate(float rate) {_rate=rate;}
    
    /// Sets the RNode pointer associated with this ReactionBase to rhs. Usually is called only by the
    /// Gillespie-like algorithms.
    void setRnode(RNode *rhs) {_rnode=rhs;}
    
    /// Returns the rate associated with this ReactionBase.
    float getRate() const {return _rate;}
    
    /// Returns the bare rate associated with this ReacitonBase
    float getBareRate() const {return _rate_bare;}
    
    ///For cross filament reactions (SHOULD EVENTUALLY BE MOVED)
    
    /// Set rMin for linker and motor reactions
    void setRMin(double rMin) {_rMin = rMin;}
    /// Set rMax for linker and motor reactions
    void setRMax(double rMax) {_rMax = rMax;}
    
    ///getters for rMin and rMax
    double getRMin() {return _rMin;}
    double getRMax() {return _rMax;}
    
    ///Set and get ID
    void setReactionID(int ID) {_reactionID = ID;}
    int getReactionID() {return _reactionID;}
    
    /// Returns a pointer to the RNode associated with this ReactionBase.
    RNode* getRnode() const {return _rnode;}
    
    /// Returns the number of reactant RSpecies
    unsigned short getM() const {return getMImpl();}
    
    /// (Private) implementation of the getM() method to be elaborated in derived classes.
    virtual unsigned short getMImpl() const = 0;
    
    /// Returns the number of product RSpecies
    unsigned short getN() const {return getNImpl();}
    
    /// (Private) implementation of the getN() method to be elaborated in derived classes.
    virtual unsigned short getNImpl() const = 0;
    
    /// Returns the total number of RSpecies
    unsigned short size() const {return sizeImpl();}
    
    /// (Private) implementation of the size() method to be elaborated in derived classes.
    virtual unsigned short sizeImpl() const = 0;
    
    /// Return the parent Composite object pointer, to which this Reaction belongs to. If not present,
    /// return a nullptr.
    Composite* getParent() {return _parent;}
    
    /// Set the parent Composite object pointer to which this Reaction belongs to.
    void setParent (Composite *other) {_parent=other;}
    
    /// Returns true if this Reaction has a parent Composite object to which it belongs to.
    bool hasParent() const {return _parent!=nullptr? true : false;}
    
    /// Get the root parent (i.e. follow the pointers of parentage until the root node in the Composition hieararchy)
    Composite* getRoot();
    
    //    Composite* getParent() {return nullptr;}
    //
    //    void setParent (Composite *other) {}
    //
    //    bool hasParent() const {return false;}
    //
    //    Composite* getRoot() {return nullptr;}
    
    /// Computes the product of the copy number of all reactant RSpecies.
    /// Can be used to quickly determine whether this ReactionBase should be allowed to activate - if one of the
    /// reactants has a copy number equal to zero, then zero is returned, indicating that
    /// this ReactionBase should not be (yet) activated.
    int getProductOfReactants () const {return getProductOfReactantsImpl();}
    
    /// (Private) implementation of the getProductOfReactants() method to be elaborated in derived classes.
    virtual int getProductOfReactantsImpl() const = 0;
    
    /// Computes the product of the copy number of all product RSpecies minus maxium allowed copy number.
    /// Can be used to quickly determine whether this ReactionBase should be allowed to activate - if one of the
    /// products has a copy number equal to the maximum allowed, then zero is returned, indicating that
    /// this ReactionBase should not be (yet) activated.
    int getProductOfProducts () const {return getProductOfProductsImpl();}
    
    /// (Private) implementation of the getProductOfProducts() method to be elaborated in derived classes.
    virtual int getProductOfProductsImpl() const = 0;
    
    /// Return true if the ReactionBase is currently passivated
#if defined TRACK_ZERO_COPY_N || defined  TRACK_UPPER_COPY_N
    bool isPassivated() const {return _passivated;}
#else
    bool isPassivated() const {return false;}
#endif
    
    /// Returns true of this Reaction contains Species *s either as a reactant or a product
    bool containsSpecies(Species *s) const {return containsSpeciesImpl(s);}
    
    /// (Private) implementation of the containsSpecies() method to be elaborated in derived classes.
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
    
#ifdef RSPECIES_SIGNALING
    virtual void broadcastRSpeciesSignals() = 0;
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
    
    /// (Private) implementation of the operator==(...) method to be elaborated in derived classes.
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
    
    /// (Private) implementation of the makeStep() method to be elaborated in derived classes.
    virtual void makeStepImpl() = 0;
    
    /// Compute the ReactionBase propensity that is needed by a Gillespie like algorithm:
    /// rate*reactant_1.getN()*reactant_2.getN()...
    float computePropensity() const {return computePropensityImpl();}
    
    /// (Private) implementation of the computePropensity() method to be elaborated in derived classes.
    virtual float computePropensityImpl() const = 0;
    
    /// Usually is applied to ReactionBase objects with propensity of 0 (e.g. when one of the
    /// copy numbers of reactants has dropped to 0. This method call notifies all other
    /// ReactionBase objects that may affect this ReactionBase to stop tracking this ReactionBase. Eventually,
    /// activateReaction() may be called to restart tracking, if the propensity stops being 0.
    void passivateReaction() {passivateReactionImpl();}
    
    /// (Private) implementation of the passivateReaction() method to be elaborated in derived classes.
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
    
//        
//        ///Replace a given species with a new one. Registers and unregisters new dependents accordingly
//        virtual void replaceRSpecies(RSpecies* oldRSpecies, RSpecies* newRSpecies) = 0;
    
    virtual void printToStream(std::ostream& os) const = 0;
    
    /// Print self into an iostream
    friend std::ostream& operator<<(std::ostream& os, const ReactionBase& rr)
    {
        rr.printToStream(os);
        return os;
    }
};



#endif /* defined(__CytoSim__ReactionBase__) */
