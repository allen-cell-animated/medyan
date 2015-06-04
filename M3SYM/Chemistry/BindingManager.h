
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

#ifndef M3SYM_BindingManager_h
#define M3SYM_BindingManager_h

#include <unordered_set>
#include <unordered_map>
#include <random>

#include "common.h"

#include "NeighborListImpl.h"
#include "ReactionBase.h"

#include "SysParams.h"

#define SPECIESBB_BINDING_INDEX 2
#define SPECIESMLB_BINDING_INDEX 0

//FORWARD DECLARATIONS
class ReactionBase;
class CCylinder;
class Compartment;

/// To store and manage binding reactions.

/*!
 *  FilamentBindingManager is used to store a binding Reaction on [Filaments](@ref Filament) 
 *  in Compartment. Contains the binding reaction, possible binding sites, and integers 
 *  representing the binding species involved in the reaction. Classes that extend this will 
 *  implement their own data structures for holding possible reaction sites, etc.
 *
 *  The main function of this class is to call updatePossibleBindings(), which will
 *  update the possible binding sites if the binding reaction is called in this compartment.
 *  Also contains functions to replace ccylinders in the structure, as well as remove a ccylinder
 *  from the structure when it is removed from the subsystem.
 *
 *  When the binding reaction is fired, the manager will choose a random binding based on its current
 *  data structure state, and perform whatever callback is associated with that binding reaction.
 */
class FilamentBindingManager {
    
friend class ChemManager;
    
protected:
    ReactionBase* _bindingReaction; ///< The binding reaction for this compartment
    
    Compartment* _compartment; ///< Compartment this is in
    
    short _boundInt; ///< Integer index in CMonomer of bound chemical value
    string _boundName; ///< String name of bound chemical value
    
    Species* _bindingSpecies; ///< The binding species that this manager tracks.
                              ///< Resposible for all copy number changes
    
    short _nlIndex = 0; ///<Index of this manager (for access of neighbor lists)
    short _mIndex = 0;  ///<Index of this manager (for access in other compartments)
    
    static mt19937 *_eng; ///< Random number generator
    
    ///helper function to update copy number and reactions
    void updateBindingReaction(int oldN, int newN) {
        
        int diff = newN - oldN;
        
        //update copy number
        if(diff > 0) {
            while (diff != 0) {
                _bindingSpecies->up();
                diff--;
            }
        }
        else if(diff < 0) {
            while (diff != 0) {
                _bindingSpecies->down();
                diff++;
            }
        }
        else {} //do nothing
        
        //check if matching
        assert((_bindingSpecies->getN() == numBindingSites())
               && "Species representing binding sites does not match \
                   number of binding sites held by the manager.");
        
        _bindingReaction->updatePropensity();
    }
    
public:
    FilamentBindingManager(ReactionBase* reaction,
                           Compartment* compartment,
                           short boundInt, string boundName)
    
    : _bindingReaction(reaction), _compartment(compartment),
      _boundInt(boundInt), _boundName(boundName) {
    
#if !defined(REACTION_SIGNALING) || !defined(RSPECIES_SIGNALING)
        cout << "Any filament binding reaction relies on reaction and species signaling. Please"
        << " set these compilation macros and try again. Exiting." << endl;
        exit(EXIT_FAILURE);
#endif

    }
    ~FilamentBindingManager() {delete _eng;}
    
    //@{
    ///add possible binding reactions that could occur
    virtual void addPossibleBindings(CCylinder* cc, short bindingSite) = 0;
    virtual void addPossibleBindings(CCylinder* cc) = 0;
    //@}
    
    //@{
    /// Remove all bindings including this cylinder
    virtual void removePossibleBindings(CCylinder* cc, short bindingSite) = 0;
    virtual void removePossibleBindings(CCylinder* cc) = 0;
    //@}
    
    ///update all possible binding reactions that could occur
    virtual void updateAllPossibleBindings() = 0;
    
    /// Get current number of binding sites
    virtual int numBindingSites() = 0;
    
    ///Get the bound species integer index
    short getBoundInt() {return _boundInt;}
    ///Get the bound species name
    string getBoundName() {return _boundName;}
    
    ///Set the index of this manager, for access to NeighborList
    void setNLIndex(int index) {_nlIndex = index;}
    
    ///Set the index of this manager, for access to other managers
    void setMIndex(int index) {_mIndex = index;}
};


/// Manager for Filament and BranchPoint creation
class BranchingManager : public FilamentBindingManager {

friend class ChemManager;
    
private:
    ///possible bindings at current state
    unordered_set<tuple<CCylinder*, short>> _possibleBindings;
    
public:
    BranchingManager(ReactionBase* reaction,
                     Compartment* compartment,
                     short boundInt, string boundName);
    ~BranchingManager() {}
    
    //@{
    ///add possible binding reactions that could occur
    virtual void addPossibleBindings(CCylinder* cc, short bindingSite);
    virtual void addPossibleBindings(CCylinder* cc);
    //@}
    
    //@{
    /// Remove all bindings including this cylinder
    virtual void removePossibleBindings(CCylinder* cc, short bindingSite);
    virtual void removePossibleBindings(CCylinder* cc);
    //@}
    
    ///update all possible binding reactions that could occur
    virtual void updateAllPossibleBindings();
    
    virtual int numBindingSites() {
        
        return _possibleBindings.size();
    }
    
    /// Choose a random binding site based on current state
    tuple<CCylinder*, short> chooseBindingSite() {
        
        assert((_possibleBindings.size() != 0)
               && "Major bug: Branching manager should not have zero binding \
                  sites when called to choose a binding site.");
        
        uniform_int_distribution<> dis(0, _possibleBindings.size() - 1);
        
        int randomIndex = dis(*_eng);
        auto it = _possibleBindings.begin();
        
        advance(it, randomIndex);
        
        return *it;
    }
};

/// Manager for Linker binding.
/*!
 *  LinkerBindingManager controls linker binding in a compartment.
 *  Manages a multimap of possible binding sites.
 */
class LinkerBindingManager : public FilamentBindingManager {
    
friend class ChemManager;
    
private:
    float _rMin; ///< Minimum reaction range
    float _rMax; ///< Maximum reaction range
    
    //possible bindings at current state. updated according to neighbor list
    unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>> _possibleBindings;
    
    //static neighbor list
    static vector<CCNeighborList*> _neighborLists;
    
public:
    LinkerBindingManager(ReactionBase* reaction,
                         Compartment* compartment,
                         short boundInt, string boundName,
                         float rMax, float rMin);
    
    ~LinkerBindingManager() {}
    
    //@{
    ///add possible binding reactions that could occur
    virtual void addPossibleBindings(CCylinder* cc, short bindingSite);
    virtual void addPossibleBindings(CCylinder* cc);
    //@}
    
    //@{
    /// Remove all bindings including this cylinder
    virtual void removePossibleBindings(CCylinder* cc, short bindingSite);
    virtual void removePossibleBindings(CCylinder* cc);
    //@}
    
    ///update all possible binding reactions that could occur
    virtual void updateAllPossibleBindings();
    
    virtual int numBindingSites() {
        
        return _possibleBindings.size();
    }
    
    //@{
    /// Getters for distances
    float getRMin() {return _rMin;}
    float getRMax() {return _rMax;}
    //@}
    
    /// Choose random binding sites based on current state
    vector<tuple<CCylinder*, short>> chooseBindingSites() {
        
        assert((_possibleBindings.size() != 0)
               && "Major bug: Linker binding manager should not have zero binding \
                   sites when called to choose a binding site.");
        
        uniform_int_distribution<> dis(0, _possibleBindings.size() - 1);
        
        int randomIndex = dis(*_eng);
        auto it = _possibleBindings.begin();
        
        advance(it, randomIndex);
        
        return vector<tuple<CCylinder*, short>>{it->first, it->second};
    }
};

/// Manager for MotorGhost binding
/*!
 *  MotorBindingManager controls motor binding in a compartment.
 *  Manages a multimap of possible binding sites.
 */
class MotorBindingManager : public FilamentBindingManager {
    
friend class ChemManager;
    
private:
    float _rMin; ///< Minimum reaction range
    float _rMax; ///< Maximum reaction range
    
    //possible bindings at current state. updated according to neighbor list
    unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>> _possibleBindings;
    
    //static neighbor list
    static vector<CCNeighborList*> _neighborLists;
    
public:
    MotorBindingManager(ReactionBase* reaction,
                        Compartment* compartment,
                        short boundInt, string boundName,
                        float rMax, float rMin);
    
    ~MotorBindingManager() {}
    
    //@{
    ///add possible binding reactions that could occur
    virtual void addPossibleBindings(CCylinder* cc, short bindingSite);
    virtual void addPossibleBindings(CCylinder* cc);
    //@}
    
    //@{
    /// Remove all bindings including this cylinder
    virtual void removePossibleBindings(CCylinder* cc, short bindingSite);
    virtual void removePossibleBindings(CCylinder* cc);
    //@}
    
    ///update all possible binding reactions that could occur
    virtual void updateAllPossibleBindings();
    
    virtual int numBindingSites() {
        
        return _possibleBindings.size();
    }
    
    //@{
    /// Getters for distances
    float getRMin() {return _rMin;}
    float getRMax() {return _rMax;}
    //@}
    
    /// Choose random binding sites based on current state
    vector<tuple<CCylinder*, short>> chooseBindingSites() {
        
        assert((_possibleBindings.size() != 0)
               && "Major bug: Motor binding manager should not have zero binding \
                   sites when called to choose a binding site.");
        
        uniform_int_distribution<> dis(0, _possibleBindings.size() - 1);
        
        int randomIndex = dis(*_eng);
        auto it = _possibleBindings.begin();
        
        advance(it, randomIndex);
        
        return vector<tuple<CCylinder*, short>>{it->first, it->second};
    }
};

#endif
