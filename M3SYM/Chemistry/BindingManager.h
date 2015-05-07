
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

#include "NeighborListContainer.h"
#include "ReactionBase.h"

#include "SysParams.h"

//FORWARD DECLARATIONS
class ReactionBase;
class CCylinder;
class Compartment;

/// To store and manager binding reactions.

/*!
 *  FilamentBindingManager is used to store a binding reaction on filaments in compartments.
 *  Contains the binding reaction, possible binding sites, and integers representing the
 *  binding species involved in the reaction. Classes that extend this will implement their
 *  own data structures for holding possible reaction sites, etc.
 *
 *  The main function of this class is to call updatePossibleBindings(), which will
 *  update the possible binding sites if the binding reaction is called in this compartment.
 *  Also contains functions to replace ccylinders in the structure, as well as remove a ccylinder
 *  from the structure when it is removed from the subsystem.
 *
 *  When the reaction is fired, the manager will choose a random binding based on its current
 *  data structure state, and perform whatever callback is associated with that binding reaction.
 */
class FilamentBindingManager {
    
    friend class SimpleManagerImpl;
    
protected:
    ReactionBase* _bindingReaction; ///< The binding reaction for this compartment
    
    Compartment* _compartment; ///< Compartment this is in
    
    short _boundInt; ///< Integer index of bound chemical value
    string _boundName; ///< String name of bound chemical value
    
    short _index = 0; ///<Index of this manager (for access of neighbor lists)
    
    static mt19937 *_eng; ///< Random number generator
    
public:
    FilamentBindingManager(ReactionBase* reaction, Compartment* compartment,
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
    
    ///update the possible binding reactions that could occur
    virtual void updatePossibleBindings(CCylinder* cc, short bindingSite) = 0;
    
    /// In the case of a cylinder copy, copy all possible bindings
    /// to a new cylinder.
    virtual void replacePossibleBindings(CCylinder* oldcc, CCylinder* newcc) = 0;
    
    ///update all possible binding reactions that could occur
    virtual void updateAllPossibleBindings() = 0;
    
    /// Get current number of binding sites
    virtual int numBindingSites() = 0;
    
    //update the binding reaction (called due to a copy number change)
    void updateBindingReaction() {_bindingReaction->activateReaction();}
    
    ///Get the bound species integer index
    short getBoundInt() {return _boundInt;}
    ///Get the bound species name
    string getBoundName() {return _boundName;}
    
    ///Set the index of this manager
    void setIndex(int index) {_index = index;}
};


/// Manager for Filament and BranchPoint creation
class BranchingManager : public FilamentBindingManager {
    
private:
    ///possible bindings at current state
    unordered_set<tuple<CCylinder*, short>> _possibleBindings;
    
public:
    BranchingManager(ReactionBase* reaction, Compartment* compartment,
                     short boundInt, string boundName)
    : FilamentBindingManager(reaction, compartment, boundInt, boundName) {}
    
    ~BranchingManager() {}
    
    virtual void updatePossibleBindings(CCylinder* cc, short bindingSite);
    
    virtual void replacePossibleBindings(CCylinder* oldcc, CCylinder* newcc);
    
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
 *  Also a subclass of CCNLContainer, contains a cylinder neighbors list of
 *  cylinders within range of binding for this reaction
 */
class LinkerBindingManager : public FilamentBindingManager {
    
friend class SimpleManagerImpl;
    
private:
    float _rMin; ///< Minimum reaction range
    float _rMax; ///< Maximum reaction range
    
    //possible bindings at current state. updated according to neighbor list
    unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>> _possibleBindings;
    
    //static neighbor list
    static vector<CCNLContainer*> _nlContainers;
    
public:
    LinkerBindingManager(ReactionBase* reaction, Compartment* compartment,
                         short boundInt, string boundName, float rMax, float rMin)
    
    : FilamentBindingManager(reaction, compartment, boundInt, boundName),
      _rMin(rMin), _rMax(rMax) {}
    
    ~LinkerBindingManager() {}
    
    virtual void updatePossibleBindings(CCylinder* cc, short bindingSite);
    
    virtual void replacePossibleBindings(CCylinder* oldcc, CCylinder* newcc);
    
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
 *  Also a subclass of CCNLContainer, contains a cylinder neighbors list of
 *  cylinders within range of binding for this reaction
 */
class MotorBindingManager : public FilamentBindingManager {
    
friend class SimpleManagerImpl;
    
private:
    float _rMin; ///< Minimum reaction range
    float _rMax; ///< Maximum reaction range
    
    //possible bindings at current state. updated according to neighbor list
    unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>> _possibleBindings;
    
    //static neighbor list
    static vector<CCNLContainer*> _nlContainers;
    
public:
    MotorBindingManager(ReactionBase* reaction, Compartment* compartment,
                        short boundInt, string boundName, float rMax, float rMin)
    
    : FilamentBindingManager(reaction, compartment, boundInt, boundName),
      _rMin(rMin), _rMax(rMax) {}
    
    ~MotorBindingManager() {}
    
    virtual void updatePossibleBindings(CCylinder* cc, short bindingSite);
    
    virtual void replacePossibleBindings(CCylinder* oldcc, CCylinder* newcc);
    
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
