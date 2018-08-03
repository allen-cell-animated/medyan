
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_BindingManager_h
#define MEDYAN_BindingManager_h

#include <unordered_set>
#include <unordered_map>
#include <random>

#include "common.h"

#include "NeighborListImpl.h"
#include "ReactionBase.h"

#include "SysParams.h"
#include "Rand.h"

#define B_RXN_INDEX 2
#define ML_RXN_INDEX 0

//FORWARD DECLARATIONS
class SubSystem;
class ReactionBase;
class CCylinder;
class Compartment;
class Cylinder;

///Enumeration for nucleation zone type. Used by CaMKIIingManager.
enum NucleationZoneType {
    ALL, BOUNDARY, TOPBOUNDARY
};

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
    
    short _boundInt; ///< Integer index in CMonomer of bound chemical value.
                     ///< @note - THIS ALSO REPRESENTS THE SPECIES TYPE THAT IS MANAGED.
    
    string _boundName; ///< String name of bound chemical value
    
    short _filamentType; ///< The filament type to operate on
    
    Species* _bindingSpecies; ///< The binding species that this manager tracks.
                              ///< Resposible for all copy number changes
    
    short _nlIndex = 0; ///<Index of this manager (for access of neighbor lists)
    short _mIndex = 0;  ///<Index of this manager (for access in other compartments)
    
    static SubSystem *_subSystem; ///< Ptr to the SubSystem
    
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
                           short boundInt, string boundName,
                           short filamentType)
    
    : _bindingReaction(reaction), _compartment(compartment),
      _boundInt(boundInt), _boundName(boundName), _filamentType(filamentType) {
    
#if !defined(REACTION_SIGNALING) || !defined(RSPECIES_SIGNALING)
        cout << "Any filament binding reaction relies on reaction and species signaling. Please"
        << " set these compilation macros and try again. Exiting." << endl;
        exit(EXIT_FAILURE);
#endif

    }
    ~FilamentBindingManager() {}
    
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
    
    ///Check consistency and correctness of binding sites. Used for debugging.
    virtual bool isConsistent() = 0;
    
    ///get the filament that the species binds to aravind June 30, 2016.
    short getfilamentType() {return _filamentType;}
    
    /// ARAVIND ADDED FEB 17 2016. append possible bindings.
    virtual void appendpossibleBindings(tuple<CCylinder*, short> t1, tuple<CCylinder*, short> t2)=0;
    ///aravind, June 30,2016.
    vector<string> getrxnspecies(){return _bindingReaction->getreactantspecies();}
    virtual void clearpossibleBindings()=0;
};


/// Manager for Filament and CaMKIIingPoint creation
class CaMKIIingManager : public FilamentBindingManager {

friend class ChemManager;
    
private:
    ///Nucleation zone type, to define where nucleation should occur
    NucleationZoneType _nucleationZone;
    
    ///If using a nucleation zone, nucleating distance from the boundary
    double _nucleationDistance;
    
    ///possible bindings at current state
    unordered_set<tuple<CCylinder*, short>> _possibleBindings;
    vector<tuple<tuple<CCylinder*, short>, tuple<CCylinder*, short>>> _camkiirestarttuple; //Used only during restart conditions.
public:
    CaMKIIingManager(ReactionBase* reaction,
                     Compartment* compartment,
                     short boundInt, string boundName,
                     short filamentType,
                     NucleationZoneType zone = NucleationZoneType::ALL,
                     double nucleationDistance = numeric_limits<double>::infinity());
    ~CaMKIIingManager() {}
    
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
               && "Major bug: CaMKIIing manager should not have zero binding \
                  sites when called to choose a binding site.");
        
        int randomIndex = Rand::randInteger(0, _possibleBindings.size() - 1);
        auto it = _possibleBindings.begin();
        
        advance(it, randomIndex);
        
        return *it;
    }
    
    virtual bool isConsistent();
    /// ARAVIND ADDED FEB 17 2016. append possible bindings.
    virtual void appendpossibleBindings(tuple<CCylinder*, short> t1, tuple<CCylinder*, short> t2){
        double oldN=numBindingSites();
        _possibleBindings.insert(t1);
        _camkiirestarttuple.push_back(make_tuple(t1,t2));
//        _camkiiCylinder=(get<0>(t2));
        double newN=numBindingSites();
        updateBindingReaction(oldN,newN);}
    virtual void clearpossibleBindings() {
        double oldN=numBindingSites();
        _possibleBindings.clear();
        updateBindingReaction(oldN,0);
    }
    vector<tuple<tuple<CCylinder*, short>, tuple<CCylinder*, short>>> getbtuple() {
        return _camkiirestarttuple;
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
    static vector<CylinderCylinderNL*> _neighborLists;
    
public:
    LinkerBindingManager(ReactionBase* reaction,
                         Compartment* compartment,
                         short boundInt,
                         string boundName,
                         short filamentType,
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
    /// ARAVIND ADDED FEB 17 2016. append possible bindings.
    virtual void appendpossibleBindings(tuple<CCylinder*, short> t1, tuple<CCylinder*, short> t2){
        double oldN=numBindingSites();
         _possibleBindings.emplace(t1,t2);
         double newN=numBindingSites();
        updateBindingReaction(oldN,newN);
    }
    //get possible bindings.
    virtual unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>> getpossibleBindings(){
        return _possibleBindings;
    }
    
    /// ARAVIND ADDED FEB 18 2016. clear possible bindings.
    virtual void clearpossibleBindings() {
        double oldN=numBindingSites();
        _possibleBindings.clear();
        updateBindingReaction(oldN,0);
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
        
        int randomIndex = Rand::randInteger(0, _possibleBindings.size() - 1);
        auto it = _possibleBindings.begin();
        
        advance(it, randomIndex);
        
        return vector<tuple<CCylinder*, short>>{it->first, it->second};
    }
    
    virtual bool isConsistent();
};

/// Manager for MotorGhost binding
/*!
 *  MotorBindingManager controls motor binding in a compartment.
 *  Manages a multimap of possible binding sites.
 */
class MotorBindingManager : public FilamentBindingManager {
    
friend class ChemManager;
    
private:
        //DEPRECATED AS OF 9/22/16
//    vector<int> _unboundIDs = {};
    ///< A vector of unbound motor ID's that are contained in this compartment. This is used
    ///< for tracking binding/unbinding and movement of specific motors.
    
    float _rMin; ///< Minimum reaction range
    float _rMax; ///< Maximum reaction range
    
    //possible bindings at current state. updated according to neighbor list
    unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>> _possibleBindings;
    
    //static neighbor list
    static vector<CylinderCylinderNL*> _neighborLists;
    
public:
    MotorBindingManager(ReactionBase* reaction,
                        Compartment* compartment,
                        short boundInt,
                        string boundName,
                        short filamentType,
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
    /// ARAVIND ADDED FEB 22 2016. append possible bindings.
    virtual void appendpossibleBindings(tuple<CCylinder*, short> t1, tuple<CCylinder*, short> t2){
        double oldN=numBindingSites();
        _possibleBindings.emplace(t1,t2);
        double newN=numBindingSites();
        updateBindingReaction(oldN,newN);
    }
    //get possible bindings.
    virtual unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>> getpossibleBindings(){
        return _possibleBindings;
    }
    
    /// ARAVIND ADDED FEB 18 2016. clear possible bindings.
    virtual void clearpossibleBindings() {
        double oldN=numBindingSites();
        _possibleBindings.clear();
        updateBindingReaction(oldN,0);
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
        
        int randomIndex = Rand::randInteger(0, _possibleBindings.size() - 1);
        auto it = _possibleBindings.begin();
        
        advance(it, randomIndex);
        
        return vector<tuple<CCylinder*, short>>{it->first, it->second};
    }
    
    virtual bool isConsistent();
    
    
    //DEPRECATED AS OF 9/8/16
//    /// Adds an unbound ID to the container
//    void addUnboundID(int ID) {_unboundIDs.push_back(ID);}
//    
//    /// Get a random ID from the container, and remove the ID
//    int getUnboundID() {
//        
//        assert(_unboundIDs.size() != 0 && "Major bug: No unbound IDs, but non-zero copy numbers.");
//        
//        int ri = Rand::randInteger(0, _unboundIDs.size() - 1);
//        int ID = _unboundIDs[ri];
//        
//        //delete and return
//        _unboundIDs.erase(_unboundIDs.begin() + ri);
//        
//        return ID;
//    }
//    
//    // remove a specific ID from the list
//    void removeUnboundID(int ID) {
//        
//        for (auto IDit = _unboundIDs.begin(); IDit != _unboundIDs.end(); IDit++) {
//            
//            if(*IDit == ID) {
//                _unboundIDs.erase(IDit);
//                return;
//            }
//        }
//    }
//    
//    ///Get all unbound ID's, but do not change container
//    const vector<int>& getAllUnboundIDs() const { return _unboundIDs; }

};

#endif