
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

#ifndef M3SYM__FilamentBindingManager_h
#define M3SYM__FilamentBindingManager_h

#include <unordered_set>
#include <unordered_map>
#include <random>

#include "common.h"

#include "NeighborListContainer.h"

#include "SysParams.h"

//FORWARD DECLARATIONS
class ReactionBase;
class CCylinder;

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
    
    /// Integer index of bound chemical value
    short _bound;
    
    static mt19937 *_eng; ///< Random number generator
    
public:
    FilamentBindingManager(ReactionBase* reaction, short bound)
    
    : _bindingReaction(reaction), _bound(bound) {
    
#if !defined(REACTION_SIGNALING)
        cout << "Any filament binding reaction relies on reaction signaling. Please"
        << " set this compilation macro and try again. Exiting." << endl;
        exit(EXIT_FAILURE);
#endif
    }
    ~FilamentBindingManager() {delete _eng;}
    
    ///update the possible binding reactions that could occur
    virtual void updatePossibleBindings(CCylinder* cc) = 0;
    
    ///remove possible binding reactions for a given ccylinder
    virtual void removePossibleBindings(CCylinder* cc) = 0;
    
    /// In the case of a cylinder copy, copy all possible bindings
    /// to a new cylinder.
    virtual void replacePossibleBindings(CCylinder* oldcc, CCylinder* newcc) = 0;
    
    ///Get the bound species
    virtual short getBound() {return _bound;}
};


/// Manager for Filament and BranchPoint creation
class BranchingManager : public FilamentBindingManager {
    
private:
    ///possible bindings at current state
    unordered_set<tuple<CCylinder*, short>> _possibleBindings;
    
public:
    BranchingManager(ReactionBase* reaction, short bound)
    : FilamentBindingManager(reaction, bound) {}
    
    ~BranchingManager() {}
    
    virtual void updatePossibleBindings(CCylinder* cc);
    virtual void removePossibleBindings(CCylinder* cc);
    virtual void replacePossibleBindings(CCylinder* oldcc, CCylinder* newcc);
    
    /// Choose a random binding site based on current state
    tuple<CCylinder*, short> chooseBindingSite() {
        
        uniform_int_distribution<> dis(0, _possibleBindings.size() - 1);
        
        int randomIndex = dis(_eng);
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
class LinkerBindingManager : public FilamentBindingManager, public CCNLContainer {
    
private:
    
    float _rMin; ///< Minimum reaction range
    float _rMax; ///< Maximum reaction range
    
    //possible bindings at current state. updated according to neighbor list
    unordered_map<tuple<CCylinder*, short>, tuple<CCylinder*, short>> _possibleBindings;
    
public:
    LinkerBindingManager(ReactionBase* reaction, short bound, float rMax, float rMin)
    
    : FilamentBindingManager(reaction, bound),
    CCNLContainer(rMax + SysParams::Geometry().cylinderSize,
                  max(rMin - SysParams::Geometry().cylinderSize, 0.0)),
    
    _rMin(rMin), _rMax(rMax) {}
    
    ~LinkerBindingManager() {}
    
    virtual void updatePossibleBindings(CCylinder* cc);
    virtual void removePossibleBindings(CCylinder* cc);
    virtual void replacePossibleBindings(CCylinder* oldcc, CCylinder* newcc);
    
    /// Choose random binding sites based on current state
    vector<tuple<CCylinder*, short>> chooseBindingSites() {
        
        uniform_int_distribution<> dis(0, _possibleBindings.size() - 1);
        
        int randomIndex = dis(_eng);
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
class MotorBindingManager : public FilamentBindingManager, public CCNLContainer {
    
private:
    float _rMin; ///< Minimum reaction range
    float _rMax; ///< Maximum reaction range
    
    //possible bindings at current state. updated according to neighbor list
    unordered_map<tuple<CCylinder*, short>, tuple<CCylinder*, short>> _possibleBindings;
    
public:
    MotorBindingManager(ReactionBase* reaction, short bound, float rMax, float rMin)
    
    : FilamentBindingManager(reaction, bound),
    CCNLContainer(rMax + SysParams::Geometry().cylinderSize,
                  max(rMin - SysParams::Geometry().cylinderSize, 0.0)),
    
    _rMin(rMin), _rMax(rMax) {}
    
    ~MotorBindingManager() {}
    
    virtual void updatePossibleBindings(CCylinder* cc);
    virtual void removePossibleBindings(CCylinder* cc);
    virtual void replacePossibleBindings(CCylinder* oldcc, CCylinder* newcc);
    
    /// Choose random binding sites based on current state
    vector<tuple<CCylinder*, short>> chooseBindingSites() {
        
        uniform_int_distribution<> dis(0, _possibleBindings.size() - 1);
        
        int randomIndex = dis(_eng);
        auto it = _possibleBindings.begin();
        
        advance(it, randomIndex);
        
        return vector<tuple<CCylinder*, short>>{it->first, it->second};
    }
};

#endif
