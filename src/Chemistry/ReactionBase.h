
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_ReactionBase_h
#define MEDYAN_ReactionBase_h

#include <iostream>

#include <iostream>
#include <array>
#include <algorithm>
#include <utility>
#include <set>
#include <unordered_set>
#include <cmath>
#include <initializer_list>

#include <boost/signals2/signal.hpp>
#include <boost/signals2/connection.hpp>
#include <boost/signals2/shared_connection_block.hpp>

#include "common.h"

#include "Species.h"

//FORWARD DECLARATIONS
class CBound;
class RNode;
class Composite;
class ReactionBase;
class SpeciesPtrContainerVector;

///Enumeration for type of reaction
enum ReactionType {
    REGULAR, DIFFUSION, POLYMERIZATIONPLUSEND, POLYMERIZATIONMINUSEND,
    DEPOLYMERIZATIONPLUSEND, DEPOLYMERIZATIONMINUSEND,
    LINKERBINDING, MOTORBINDING, LINKERUNBINDING, MOTORUNBINDING,
    MOTORWALKINGFORWARD, MOTORWALKINGBACKWARD,
    AGING, FILAMENTCREATION, FILAMENTDESTRUCTION,
    BRANCHING, BRANCHUNBINDING, SEVERING
};

/// This is a ReactionBase signal object that may be called by a ReactionBase
/// simulation algorithm
typedef boost::signals2::signal<void (ReactionBase *)> ReactionEventSignal;


/// Represents an abstract interface for simple chemical reactions of the form
/// A + B -> C.

/*! ReactionBase provides an interface for managing a chemical reaction. It is an 
 *  abstract interface, so it cannot be directly instantiated, but other concrete 
 *  classes may be derived from it. ReactionBase may have a composite object as a 
 *  parent. A signaling interface may be used to make callbacks when some event, such 
 *  as a single reaction step, has been executed.
 *  The ReactionBase indicates a forward process only. For processes in both directions, 
 *  e.g. A <-> B, two ReactionBases objects need to be defined, corresponding to A->B 
 *  and B->A.
 *  A ReactionBase tracks other ReactionBase objects that are affected if this 
 *  ReactionBase is executed. A ReactionBase may be set up such that it "signals" when a 
 *  ReactionBase event happens, in which case the corresponding callbacks are called.
 */
class ReactionBase {
protected:

    unordered_set<ReactionBase*>  _dependents; ///< Pointers to ReactionBase objects that depend
                                               ///< on this ReactionBase being executed
    
    RNode* _rnode; ///< A pointer to an RNode object which is used
                   ///< to implement a Gillespie-like algorithm (e.g. NRM)
    
    Composite *_parent; ///< A pointer to a Composite object to which
                        ///< this Reaction belongs
    
    float _rate;      ///< the rate for this ReactionBase
    float _rate_bare; ///< the bare rate for this ReactionBase (original rate)

#ifdef REACTION_SIGNALING
    unique_ptr<ReactionEventSignal> _signal;///< Can be used to broadcast a signal
                                            ///< associated with this ReactionBase
#endif
#if defined TRACK_ZERO_COPY_N || defined TRACK_UPPER_COPY_N
    bool _passivated; ///< Indicates whether the ReactionBase is currently passivated
#endif
    ReactionType _reactionType; ///< Reaction type enumeration
    
    bool _isProtoCompartment = false;///< Reaction is in proto compartment
    ///< (Do not copy as a dependent, not in ChemSim)
    
    CBound* _cBound = nullptr; ///< CBound that is attached to this reaction
    
    
    float _gnum = 0.0;
    
    string _hrcdid = "DNT";
    
    float _linkRateForward = 0.0;
    
    float _linkRateBackward = 0.0;

    floatingpoint _volumeFrac; ///< Used in compartments to store volume fraction of the compartment
    int _rateVolumeDepExp; ///< Exponent of rate dependency on volume
    ///< Dependence on bulk properties are NOT considered currently
public:
    
    /// The main constructor:
    /// @param rate - the rate constant (full volume) for this ReactionBase
    ReactionBase (float rate, bool isProtoCompartment, floatingpoint volumeFrac=1.0, int rateVolumeDepExp=0);
    
    /// No copying (including all derived classes)
    ReactionBase (const ReactionBase &rb) = delete;
    
    /// no assignment (including all derived classes)
    ReactionBase& operator=(ReactionBase &rb) = delete;
    
    /// Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~ReactionBase() noexcept {}
    
    /// Copy this reaction using SpeciesPtrContainerVector &spcv as a source of
    /// analogous Species.
    /// @return the cloned ReactionBase pointer.
    /// @note the receving object should take care of the memory management of the
    /// returned pointer
    ReactionBase* clone(const SpeciesPtrContainerVector &spcv) {
        return cloneImpl(spcv);
    }
    
    /// (Private) implementation of the clone() method to be elaborated in derived classes
    virtual ReactionBase* cloneImpl(const SpeciesPtrContainerVector &spcv) = 0;
    
    /// Returns a pointer to the first element of the array<RSpecies*>. This pointer can
    /// be used to iterate over RSpecies* if necessary (also by relying on getM() and
    /// size() to determine the iteration limits). The corresponding array<RSpecies*> is
    /// defined by the derived classes.
    virtual RSpecies** rspecies() = 0;
    
    //aravind, June 30, 2016.
    vector<string> getreactantspecies(){
        vector<string> returnvector;
        for(int i=0;i<2;i++){
            returnvector.push_back((*(rspecies()+i))->getSpecies().getName());
        }
        return returnvector;
    }

    
    vector<string> getReactantSpecies() {
        vector<string> returnvector;
        for(auto i=0U;i<getM();++i){
            string name = (*(rspecies()+i))->getSpecies().getName();
            string namecut = name.substr(0,name.find("-",0));
            returnvector.push_back(namecut);
        }
        return returnvector;
    }
    
    vector<string> getProductSpecies() {
        vector<string> returnvector;
        for(auto i=getM();i<size();++i){
        string name = (*(rspecies()+i))->getSpecies().getName();
        string namecut = name.substr(0,name.find("-",0));
        returnvector.push_back(namecut);
        }
        return returnvector;
    }
    
    vector<species_copy_t> getReactantCopyNumbers()  {
        vector<species_copy_t> returnvector;
        for(auto i=0U;i<getM();i++)
        {returnvector.push_back((*(rspecies()+i))->getN());}
        return returnvector;
    }
    
    vector<species_copy_t> getProductCopyNumbers()  {
        vector<species_copy_t> returnvector;
        for(auto i=getM();i<size();i++)
        {returnvector.push_back((*(rspecies()+i))->getN());}
        return returnvector;
    }

    
    
    ///Set reaction type
    void setReactionType(ReactionType rxnType) {_reactionType = rxnType;}
    
    ///Get reaction type
    ReactionType getReactionType() {return _reactionType;}
    
    void setGNumber(floatingpoint gnum) {_gnum = gnum;};

	floatingpoint getGNumber() {return _gnum;};
    
    void setHRCDID(string hrcdid) {_hrcdid = hrcdid;};
    
    string getHRCDID() {return _hrcdid;};
    
    ///Set CBound
    void setCBound(CBound* cBound) {_cBound = cBound;}
    ///Get CBound
    CBound* getCBound() {return _cBound;}
    
    /// Sets the ReactionBase rate to the parameter "rate"
    void setRate(float rate) { _rate = rate; }
    // Sets the scaled rate based on volume dependence.
    void setRateScaled(float rate) {
        // This can automatically set the "_rate" as scaled value of "rate"

        // Some possibilities of the exponent are implemented specifically to decrease the use of "pow"
        switch(_rateVolumeDepExp) {
        case 0:
            _rate = rate; break;
        case -1:
            _rate = rate / _volumeFrac; break;
        default:
            if(_volumeFrac == 1.0f) _rate = rate;
            else _rate = rate * std::pow(_volumeFrac, _rateVolumeDepExp);
            break;
        }
    }
    /// Returns the rate associated with this ReactionBase.
    float getRate() const {return _rate;}
    
    /// Returns the bare rate associated with this ReactionBase
    float getBareRate() const {return _rate_bare;}
    ///aravind June 24, 2016
    void setBareRate(float a) {_rate_bare=a;}

    /// Getter and setter for compartment volume fraction
    floatingpoint getVolumeFrac()const { return _volumeFrac; }
    void setVolumeFrac(float volumeFrac) { _volumeFrac = volumeFrac; }

    /// Sets the RNode pointer associated with this ReactionBase to rhs. Usually is
    /// called only by the Gillespie-like algorithms.
    void setRnode(RNode *rhs) {_rnode=rhs;}
    /// Get the RNode pointer
    RNode* getRNode() {return _rnode;}

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
    
    /// Return the parent Composite object pointer, to which this Reaction belongs to.
    /// If not present, return a nullptr.
    Composite* getParent() {return _parent;}
    
    /// Set the parent Composite object pointer to which this Reaction belongs to.
    void setParent (Composite *other) {_parent=other;}
    
    /// Returns true if this Reaction has a parent Composite object to which it
    /// belongs to.
    bool hasParent() const {return _parent!=nullptr? true : false;}
    
    /// Get the root parent (i.e. follow the pointers of parentage until the root node
    /// in the Composition hieararchy)
    Composite* getRoot();
    
    /// Computes the product of the copy number of all reactant RSpecies.
    /// Can be used to quickly determine whether this ReactionBase should be allowed to
    /// activate - if one of the reactants has a copy number equal to zero, then zero is
    /// returned, indicating that this ReactionBase should not be (yet) activated.
    floatingpoint getProductOfReactants () const {return getProductOfReactantsImpl();}
    
    /// (Private) implementation of the getProductOfReactants() method to be elaborated
    /// in derived classes.
    virtual floatingpoint getProductOfReactantsImpl() const = 0;
    
    /// Computes the product of the copy number of all product RSpecies minus maximum
    /// allowed copy number. Can be used to quickly determine whether this ReactionBase
    /// should be allowed to activate - if one of the products has a copy number equal
    /// to the maximum allowed, then zero is returned, indicating that this ReactionBase
    /// should not be (yet) activated.
    floatingpoint getProductOfProducts () const {return getProductOfProductsImpl();}
    
    /// (Private) implementation of the getProductOfProducts() method to be elaborated
    /// in derived classes.
    virtual floatingpoint getProductOfProductsImpl() const = 0;
    
    /// Return true if the ReactionBase is currently passivated
#if defined TRACK_ZERO_COPY_N || defined  TRACK_UPPER_COPY_N
    bool isPassivated() const {return _passivated;}
#else
    bool isPassivated() const {return false;}
#endif
    
    /// Returns true of this Reaction contains Species *s either as a reactant or a
    /// product
    bool containsSpecies(Species *s) const {return containsSpeciesImpl(s);}
    
    /// (Private) implementation of the containsSpecies() method to be elaborated in
    /// derived classes.
    virtual bool containsSpeciesImpl(Species *s) const = 0;
    
#ifdef REACTION_SIGNALING
    /// Return true if this RSpecies emits signals on copy number change
    bool isSignaling () const {return _signal!=nullptr;}
    
    /// Set the signaling behavior of this ReactionBase
    void startSignaling ();
    
    /// Destroy the signal associated with this ReactionBase; all associated slots
    /// will be destroyed @note To start signaling again, startSignaling() needs to be
    /// called
    void stopSignaling ();
    
    /// Connect the callback, react_callback to a signal corresponding to
    /// ReactionBase *r.
    /// @param function<void (ReactionBase *)> const &react_callback - a function
    /// object to be called (a slot)
    /// @param int priority - lower priority slots will be called first. Default is 5.
    /// Do not use priorities 1 and 2 unless absolutely essential.
    /// @return a connection object which can be used to later disconnect this
    /// particular slot or temporarily block it
    boost::signals2::connection
        connect(function<void (ReactionBase *)> const &react_callback, int priority=5);
    
    /// Broadcasts signal indicating that the ReactionBase event has taken place
    /// This method is only called by the code which runs the chemical dynamics (i.e.
    /// Gillespie-like algorithm)
    void emitSignal() {
        if(isSignaling())
            (*_signal)(this);
    }
#endif
    
    /// Return a const reference to the set of dependent ReactionBases
    /// @note One can obtain two different lists of affected ReactionBases:
    /// 1) via getAffectedReactionBases(), where the copy numbers do influence the
    /// dependencies, and 2) via dependents(), where dependencies stop being counted
    /// if the copy numbers of reactant species drop to 0.
    const unordered_set<ReactionBase*>& dependents() {return _dependents;}
    
    /// Returns true if two ReactionBase objects are equal.
    /// Two ReactionBase objects are equal if each of their reactants and products
    /// are equal
    friend bool operator==(const ReactionBase& a, const ReactionBase& b)
    {
        if(typeid(a) != typeid(b))
            return false;
        return a.is_equal(b);
        // Invoke virtual is_equal via derived subclass of a
        // (which should be the same as b)
    }
    
    /// (Private) implementation of the operator==(...) method to be elaborated
    /// in derived classes.
    virtual bool is_equal(const ReactionBase& b) const = 0;
    
    /// Return true if two ReactionBase are not equal.
    /// @see operator ==(const ReactionBase& a, const ReactionBase& b) above
    friend bool operator !=(const ReactionBase& a, const ReactionBase& b){
        return !(a==b);
    }
    
    /// Fire the ReactionBase - make a single step, where reactant RSpecies copy
    /// numbers are decreased by one, and the product RSpecies copy numbers are
    /// increased by one.
    /// @note This method does not send a ReactionBase event Signal. The latter is
    /// usually done from within a Gillespie-like algorithm.
    void makeStep() {makeStepImpl();}
    
    /// (Private) implementation of the makeStep() method to be elaborated in
    /// derived classes.
    virtual void makeStepImpl() = 0;
    
    /// Compute the ReactionBase propensity that is needed by a Gillespie like
    /// algorithm, rate*reactant_1.getN()*reactant_2.getN()...
    floatingpoint computePropensity() const {return computePropensityImpl();}
    
    /// (Private) implementation of the computePropensity() method to be elaborated
    /// in derived classes.
    virtual floatingpoint computePropensityImpl() const = 0;
    
    /// Usually is applied to ReactionBase objects with propensity of 0 (e.g. when one
    /// of the copy numbers of reactants has dropped to 0. This method call notifies all
    /// other ReactionBase objects that may affect this ReactionBase to stop tracking
    /// this ReactionBase. Eventually, activateReaction() may be called to restart
    /// tracking, if the propensity stops being 0.
    void passivateReaction() {passivateReactionImpl();}
    
    /// (Private) implementation of the passivateReaction() method to be elaborated in
    /// derived classes.
    virtual void passivateReactionImpl() = 0;
    
    /// Requests that ReactionBase objects that may affect this Reaction to start
    /// tracking it, which can be used to follow ReactionBase objects whose propensities
    /// change upong firing of some ReactionBase. This request is acted upon
    /// unconditionally.
    void activateReactionUnconditional() {activateReactionUnconditionalImpl();}
    
    virtual void activateReactionUnconditionalImpl() = 0;
    
    /// Requests that Reaction objects that may affect this Reaction to start tracking
    /// it, which can be used to follow Reaction objects whose propensities change upon
    /// firing of some Reaction. This request will be ignored if the Reaction's
    /// propensity is still zero.
    void activateReaction() {
#ifdef TRACK_ZERO_COPY_N
        if(areEqual(getProductOfReactants(), 0.0)) // One of the reactants is still at zero copy n,
                                                   // no need to activate yet...
            return;
#endif
#ifdef TRACK_UPPER_COPY_N
        if(areEqual(getProductOfProducts(), 0.0)) // One of the products is at the maximum allowed
                                                  //copy number, no need to activate yet...
            return;
#endif
        activateReactionUnconditional();
    }
    
    /// Performs a simple updating of the propensity of this ReactionBase. Does not change
    /// dependents, only updates the RNode if initialized.
    void updatePropensity() {updatePropensityImpl();}
    
    virtual void updatePropensityImpl() = 0;
    
    /// Print the ReactionBases that are affacted by this ReactionBase being fired
    void printDependents() ;
    
    /// Return the list of ReactionBase objects that are affected when this
    /// ReactionBase is fired
    /// @return a vector of pointers to the affected ReactionBase objects
    /// @note This method is "expensive" because it computes from scratch the
    /// dependencies. Importantly, the copy numbers of molecules do not influence the
    /// result of this function. \sa dependents()
    virtual vector<ReactionBase*> getAffectedReactions() = 0;
    
    /// Request that the ReactionBase *r adds this ReactionBase to its list of
    /// dependents which it affects.
    void registerNewDependent(ReactionBase *r);
    
    /// Request that the ReactionBase *r removes this ReactionBase from its list of
    /// dependents which it affects.
    /// This is usually requested when the ReactionBase propensity drops to zero (i.e.
    /// via passivateReactionBase()).
    void unregisterDependent(ReactionBase *r);
    
    virtual void printToStream(ostream& os) const = 0;
    
    /// Print self into an iostream
    friend ostream& operator<<(ostream& os, const ReactionBase& rr)
    {
        rr.printToStream(os);
        return os;
    }
    
    ///Whether the dependencies should be updated
    virtual bool updateDependencies() = 0;
};



#endif
