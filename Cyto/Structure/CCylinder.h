//
//  CCylinder.h
//  CytoSim
//
//  Created by James Komianos on 7/8/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __CytoSim__CCylinder__
#define __CytoSim__CCylinder__

#include <iostream>

#include "common.h"
#include "Compartment.h"
#include "CMonomer.h"
#include "SystemParameters.h"

class ChemSimReactionKey;
class Cylinder;

/// CCylinder class holds all Monomers and Bounds
/*! 
 *  The CCylinder Class is an template class that has lists of the monomers and bounds that it contains.
 *  it has functionality to print the current composition.
 *  Accessing a particular species in the CCylinder is possible as well.
 */
class CCylinder {
    
protected:
    std::vector<std::unique_ptr<CMonomer>> _monomers; ///< list of monomers in this ccylinder
    
    ///REACTION CONTAINERS
    std::vector<ReactionBase*> _reactions;///< list of reactions associated with ccylinder
    std::vector<ReactionBase*> _frontReactions; ///< list of reactions involving the next ccylinder
    std::vector<ReactionBase*> _backReactions; ///< list of reactions involving the previous ccylinder
    
    Compartment* _compartment; ///< compartment this ccylinder is in

    Cylinder* _pCylinder;
    
    short _size = SystemParameters::Geometry().cylinderSize / SystemParameters::Geometry().monomerSize;

public:
    ///Default constructor, sets compartment
    CCylinder(Compartment* c) : _compartment(c) {}
    
    /// Copy constructor
    /// @note This constructor will create a new ccylinder with different species and reactions
    /// within the compartment that is chosen as a parameter to the constructor. The copied and
    /// original ccylinders will not share reactions or species, but will be copied into a new compartment.
    CCylinder(const CCylinder& rhs, Compartment* c) : _compartment(c)
    {
        ///copy all monomers, bounds
        for(auto &m : rhs._monomers)
            _monomers.push_back(std::unique_ptr<CMonomer>(m->clone(c)));
        
        ///copy all reactions
        for(auto &r: rhs._reactions) addReaction(r->clone(c->speciesContainer()));
        for(auto &r: rhs._frontReactions) addFrontReaction(r->clone(c->speciesContainer()));
        for(auto &r: rhs._backReactions) addBackReaction(r->clone(c->speciesContainer()));
        
        ///Update and return
        this->updateReactions();
    }
    
    /// Assignment is not allowed
    CCylinder& operator=(CCylinder &rhs) = delete;
    
    ///Default destructor, explicitly removes monomers and bounds (including their species, rxns)
    ~CCylinder()
    {
        ///Remove all reactions
        for(auto &r: _reactions) removeReaction(r);
        for(auto &r: _frontReactions) removeReaction(r);
        for(auto &r: _backReactions) removeReaction(r);
        
        ///Remove all species
        for(auto &m: _monomers) {
            for(auto &s : m->speciesFilamentVector())
                _compartment->removeSpecies(s);
            for(auto &s : m->speciesBoundVector())
                _compartment->removeSpecies(s);
            for(auto &s : m->speciesPlusEndVector())
                _compartment->removeSpecies(s);
            for(auto &s : m->speciesMinusEndVector())
                _compartment->removeSpecies(s);
        }
    }
    
    ///Clone, calls copy constructor
    virtual CCylinder* clone(Compartment* c) {
        return new CCylinder(*this, c);
    }
    
    ///get filament compartment
    Compartment* getCompartment() {return _compartment;}
    
    ///set parent cylinder
    void setCylinder(Cylinder* c) {_pCylinder = c;}
    ///get parent cylinder
    Cylinder* getCylinder() {return _pCylinder;}
    
    ///Add a monomer to this CCylinder
    void addCMonomer(CMonomer* monomer) { _monomers.emplace_back(std::unique_ptr<CMonomer>(monomer));}
    ///Get monomer at an index
    ///@note no check on index
    CMonomer* getCMonomer(int index) {return _monomers[index].get();}
    
    ///Add a reaction to this CCylinder
    void addReaction(ReactionBase* r);
    ///Add a reaction with the front ccylinderto this CCylinder
    ///@note manage decides whether this reaction adding will also
    ///add to the chemsim and compartment.
    void addFrontReaction(ReactionBase* r, bool manage=false);
    ///Add a reaction with the back ccylinderto this CCylinder
    ///@note manage decides whether this reaction adding will also
    ///add to the chemsim and compartment.
    void addBackReaction(ReactionBase* r, bool manage=false);
    
    ///remove a filament reaction from this CCylinder
    ///@note no check on whether r is in the reactions list
    void removeReaction(ReactionBase* r);
    ///remove a filament reaction from this CCylinder
    ///@note no check on whether r is in the reactions list
    ///@note manage decides whether this reaction removing will also
    ///remove to the chemsim and compartment.
    void removeFrontReaction(ReactionBase* r, bool manage=false);
    ///remove a filament reaction from this CCylinder
    ///@note no check on whether r is in the reactions list
    ///@note manage decides whether this reaction removing will also
    ///remove to the chemsim and compartment.
    void removeBackReaction(ReactionBase* r, bool manage=false);
    
    ///Get list of reactions associated with this CCylinder
    std::vector<ReactionBase*>& getReactions() {return _reactions;}
    ///Get list of reactions associated with this CCylinder and next
    std::vector<ReactionBase*>& getFrontReactions() {return _frontReactions;}
    ///Get list of reactions associated with this CCylinder and previous
    std::vector<ReactionBase*>& getBackReactions() {return _backReactions;}
    
    ///Update all reactions associated with this CCylinder
    void updateReactions();
    
    ///Print CCylinder
    virtual void printCCylinder();
    
    ///get size of this ccylinder in number of monomers
    short size() {return _size;}
};



#endif /* defined(__CytoSim__CCylinder__) */
