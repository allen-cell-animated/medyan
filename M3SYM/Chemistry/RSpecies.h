
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_RSpecies_h
#define M3SYM_RSpecies_h

#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <cassert>
#include <stdexcept>

#include <boost/signals2/signal.hpp>
#include <boost/signals2/connection.hpp>
#include <boost/signals2/shared_connection_block.hpp>

#include "common.h"

//FORWARD DECLARATIONS
class Species;
class RSpecies;
class ReactionBase;
class CMonomer;

/// vr stands for vector of Reactions
typedef vector<ReactionBase*>::iterator vr_iterator; 
typedef vector<ReactionBase*>::const_iterator vr_const_iterator; 

/// vsp stands for vector of RSpecies
typedef vector<RSpecies*>::iterator vrsp_iterator; 
typedef vector<RSpecies*>::const_iterator vrsp_const_iterator; 

/// This is a RSpecies signal object that can be used to signal when the
/// copy number changes
typedef boost::signals2::signal<void (RSpecies *, int)> RSpeciesCopyNChangedSignal;

/// Represents the reactive aspect of chemical molecules. It tracks their copy
/// number and can be used in [Reactions](@ref Reaction).

/*!  This class represents the reactivity of chemical species. The name RSpecies stems
 *   from reacting species. RSpecies tracks the copy number of molecules and the 
 *   [Reactions](@ref Reaction) in which it is involed (@see Reaction).
 *   @note Each intantiation of RSpecies is unique, and hence, cannot be neither copied
 *   nor moved (C++11). This has significant implications - e.g., RSpecies cannot be 
 *   used in vector<RSpecies>. Instead, one should use either vector<RSpecies*> if not 
 *   owning the RSpecies pointers, or vector<unique_ptr<RSpecies>> if owning the 
 *   RSpecies pointers. A special allocator can be written such that dynamically 
 *   allocated RSpecies (through new), are arranged contigiously in memory.
 */
class RSpecies {
    friend Species;
    friend CMonomer;
    /// Reactions calls addAsReactant(), removeAsReactant() - which other classes
    /// should not call
private: //Variables
    vector<ReactionBase *> _as_reactants = {}; ///< a vector of [Reactions]
                                               ///< (@ref Reaction) where this RSpecies
                                               ///< is a Reactant
    vector<ReactionBase *> _as_products = {};  ///< a vector of [Reactions]
                                               ///< (@ref Reaction) where this RSpecies
                                               ///< is a Product
    Species& _species; ///< reference to the **parent** Species object
    species_copy_t _n; ///< Current copy number of this RSpecies
#ifdef TRACK_UPPER_COPY_N
    species_copy_t _ulim; ///< Upper limit for the copy number, afterwards all
                          ///< reactions leading to further accum. are turned off
#endif
#ifdef RSPECIES_SIGNALING
    RSpeciesCopyNChangedSignal *_signal; ///< Can be used to broadcast a signal
                                         ///< associated with change of n of
#endif                                   ///< this RSpecies (usually when a single step
                                         ///< of this Reaction occurs)
    
public:
    /// Constructors 
    /// @param parent - the Species object to which this RSpecies belongs
    /// @param n - copy number
    RSpecies (Species &parent, species_copy_t n=0, species_copy_t ulim=max_ulim) :
    _species(parent), _n(n)
    {
#ifdef TRACK_UPPER_COPY_N
        _ulim = ulim;
#endif
#ifdef RSPECIES_SIGNALING
        _signal=nullptr;
#endif
    }
    
    /// deleted copy constructor - each RSpecies is uniquely created by the parent
    /// Species
    RSpecies(const RSpecies &r) = delete;
    
    /// deleted move constructor - each RSpecies is uniquely created by the parent
    /// Species
    RSpecies(RSpecies &&r) = delete;
    
    /// deleted assignment operator - each RSpecies is uniquely created by the parent
    /// Species
    RSpecies& operator=(RSpecies&) = delete;
           
    /// Sets the copy number for this RSpecies. 
    /// @param n should be a non-negative number, but no checking is done in run time
    /// @note The operation does not emit any signals about the copy number change.
    void setN(species_copy_t n) {_n=n;}

#ifdef TRACK_UPPER_COPY_N        
    /// Sets the upper limit for the copy number of this RSpecies.
    /// @param ulim should be a non-negative number, but no checking is done in run time
    void setUpperLimitForN(species_copy_t ulim) {_ulim=ulim;}
    
    /// Return the upper limit for the copy number of this RSpecies
    species_copy_t getUpperLimitForN() const {return _ulim;}
#endif
    
    /// Increases the copy number by 1. If the copy number changes from 0 to 1, calls a
    /// "callback"-like method to activated previously passivated [Reactions](@ref
    /// Reaction), where this RSpecies is a Reactant.
    /// Also emits a signal with the change in copy number if attached.
    virtual inline void up() {
        _n+=1;
#ifdef TRACK_ZERO_COPY_N
        if(_n==1)
            activateAssocReactantReactions();
#endif
#ifdef TRACK_UPPER_COPY_N
        if(_n==_ulim)
            passivateAssocProductReactions();
#endif
        
#ifdef RSPECIES_SIGNALING
    if(isSignaling()) emitSignal(+1);
#endif
    }
    
    /// Decreases the copy number by 1. If the copy number changes becomes 0, calls a
    /// "callback"-like method to passivate [Reactions](@ref Reaction), where this
    /// RSpecies is a Reactant.
    /// Also emits a signal with the change in copy number if attached.
    virtual inline void down() {
#ifdef TRACK_UPPER_COPY_N
        species_copy_t prev_n = _n;
#endif
        _n-=1;
#ifdef TRACK_ZERO_COPY_N
        if(_n == 0)
            passivateAssocReactantReactions();
#endif
#ifdef TRACK_UPPER_COPY_N
        if(prev_n == _ulim)
            activateAssocProductReactions();
#endif
        
#ifdef RSPECIES_SIGNALING
    if(isSignaling()) emitSignal(-1);
#endif
    }
    
    // \internal This methods is called by the Reaction class during construction
    // of the Reaction where this RSpecies is involved as a Reactant
    void addAsReactant(ReactionBase *r){_as_reactants.push_back(r);}
    
    // \internal This methods is called by the Reaction class during construction
    // of the Reaction where this RSpecies is involved as a Product    
    void addAsProduct(ReactionBase *r){_as_products.push_back(r);}
    
    // \internal This method is called by the Reaction class during destruction
    // of the Reaction where this RSpecies is involved as a Reactant
    void removeAsReactant(const ReactionBase *r) {
        
        if(!_as_reactants.empty()) {
            auto rxit = find(_as_reactants.begin(),_as_reactants.end(),r);
            if(rxit!=_as_reactants.end()){
                _as_reactants.erase(rxit);
            }
            else {
                
                cout << "Did not find the as_reactant" << endl;
                
            }
        }
    }
    // \internal This method is called by the Reaction class during destruction
    // of the Reaction where this RSpecies is involved as a Product
    void removeAsProduct(const ReactionBase *r) {
        
        if (!_as_products.empty()) {
            auto rxit = find(_as_products.begin(),_as_products.end(),r);
            if(rxit!=_as_products.end()){
                _as_products.erase(rxit);
            }
            else {
                cout << "Did not find the as_product" << endl;
            }
        }
        
    }

    /// \internal Attempts to activate previously passivated [Reactions](@ref Reaction)
    /// where this RSpecies is involved as a Reactant. Usually, the Reaction was first
    /// passivated, for example if the RSpecies copy number of one of the reactants
    /// dropeed to zero. This attempt may not succeed if there are still other reactants
    /// in the same Reaction with zero copy count.
    void activateAssocReactantReactions();
    
    /// \internal Attempts to activate previously passivated [Reactions](@ref Reaction)
    /// where this RSpecies is involved as a Product. Usually, the Reaction was first
    /// passivated, for example if the RSpecies copy number of one of the reactants went
    /// up to max_ulim.
    void activateAssocProductReactions();
    
    /// \internal Passivates all [Reactions](@ref Reaction) where this RSpecies is among
    /// the reactants.
    void passivateAssocReactantReactions();
    
    /// \internal Passivates all [Reactions](@ref Reaction) where this RSpecies is among
    /// the products.
    void passivateAssocProductReactions();
           
#ifdef RSPECIES_SIGNALING
    /// Set the signaling behavior of this RSpecies
    void startSignaling();
    
    /// Destroy the signal associated with this RSpecies; all associated slots will be
    /// destroyed @note To start signaling again, startSignaling() needs to be called
    void stopSignaling();
#endif
    
public:
    /// It is required that all [Reactions](@ref Reaction) associated with this
    /// RSpecies are destructed before this RSpecies is destructed. Most of the time,
    /// this will occur naturally. If not, an assertion will ungracefully terminate the
    /// program.
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~RSpecies() noexcept;
    
#ifdef RSPECIES_SIGNALING
    /// Broadcasts signal indicating that the copy number of this RSpecies has changed
    /// This method should usually called by the code which runs the chemical dynamics
    /// (i.e. Gillespie-like algorithm)
    
    inline void emitSignal(int delta) {
        if(isSignaling())
            (*_signal)(this, delta);
    }
    
    /// Return true if this RSpecies emits signals on copy number change
    inline bool isSignaling() const {return _signal!=nullptr;}
#endif
    
    /// Return the current copy number of this RSpecies
    inline species_copy_t getN() const {return _n;}
    
    /// return parent Species as a reference
    inline Species& getSpecies() {return _species;}
    
    /// return parent Species as a const reference
    inline const Species& getSpecies() const {return _species;}
    
    /// Return the full name of this Species in a string format (e.g. "Arp2/3{Bulk}"
    string getFullName() const;
            
    /// Return vector<ReactionBase *>, which contains pointers to all [Reactions](@ref
    /// Reaction) where this RSpecies is involved as a Reactant
    inline vector<ReactionBase *>& reactantReactions(){return _as_reactants;}
    
    /// Return vector<ReactionBase *>, which contains pointers to all [Reactions](@ref
    /// Reaction) where this RSpecies is involved as a Product
    inline vector<ReactionBase *>& productReactions(){return _as_products;}
    
    /// Return vector<ReactionBase *>::iterator, which points to the beginning of all 
    /// [Reactions](@ref Reaction) where this RSpecies is involved as a Reactant
    inline vr_iterator beginReactantReactions() {return _as_reactants.begin();}
    
    /// Return vector<ReactionBase *>::iterator, which points to the beginning of all 
    /// [Reactions](@ref Reaction) where this RSpecies is involved as a Product
    inline vr_iterator beginProductReactions() {return _as_products.begin();}
    
    /// Return vector<ReactionBase *>::iterator, which points to the end of all 
    /// [Reactions](@ref Reaction) where this RSpecies is involved as a Reactant    
    inline vr_iterator endReactantReactions() {return _as_reactants.end();}
    
    /// Return vector<ReactionBase *>::iterator, which points to the end of all 
    /// [Reactions](@ref Reaction) where this RSpecies is involved as a Product
    inline vr_iterator endProductReactions() {return _as_products.end();}

#ifdef BOOST_MEM_POOL
    /// Advanced memory management
    void* operator new(size_t size);
    
    void operator delete(void* ptr) noexcept;
#endif
};

/// A constant RSpecies whose copy number does not change.
/*!
 *  The RSpeciesConst class represents a RSpecies that does not change copy number upon 
 *  reacting This is used for SpeciesDiffusing or SpeciesBulk species, whose copy 
 *  numbers are constant throughout the duration of the chemical simulation. Their
 *  constant qualifier is set at the Species construction, where a RSpeciesConst object 
 *  is created in place of a typical RSpecies.
 */
class RSpeciesConst : public RSpecies {
    
public:
    /// Constructors
    /// @param parent - the Species object to which this RSpeciesConst belongs
    /// @param n - copy number, will not change
    RSpeciesConst (Species &parent, species_copy_t n=0, species_copy_t ulim=max_ulim)
        : RSpecies(parent, n, ulim) {}
    /// deleted copy constructor - each RSpeciesConst is uniquely created by the parent
    /// Species
    RSpeciesConst(const RSpeciesConst &r) = delete;
    
    /// deleted move constructor - each RSpeciesConst is uniquely created by the parent
    /// Species
    RSpeciesConst(RSpeciesConst &&r) = delete;
    
    /// deleted assignment operator - each RSpeciesConst is uniquely created by the
    /// parent Species
    RSpeciesConst& operator=(RSpeciesConst&) = delete;
    
    //@{
    /// In constant species, this function does nothing. Copy numbers do not change.
    virtual inline void up() {}
    virtual inline void down() {}
    //@}
};


/// Print self into an iostream
ostream& operator<<(ostream& os, const RSpecies& s);

#endif
