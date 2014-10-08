//
//  CCylinder.cpp
//  Cyto
//
//  Created by James Komianos on 9/22/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "CCylinder.h"

CCylinder::CCylinder(const CCylinder& rhs, Compartment* c) : _compartment(c)
{
    
    ///copy all monomers, bounds
    for(auto &m : rhs._monomers)
        _monomers.push_back(std::unique_ptr<CMonomer>(m->clone(c)));
    
    ///copy all reactions
    for(auto &r: rhs._internalReactions)
        addInternalReaction(r->clone(c->speciesContainer()));
    
    for(auto it = rhs._crossCylinderReactions.begin(); it != rhs._crossCylinderReactions.end(); it++) {
        ///Copy map
        for(auto &r : it->second) addCrossCylinderReaction(it->first, r->clone(c->speciesContainer()));
    }
    ///Copy reacting cylinders
    for(auto &c : rhs._reactingCylinders) addReactingCylinder(c);
    
    

    ///Update and return
    this->updateReactions();
}

CCylinder::~CCylinder()
{
    ///Remove all reactions
    for(auto &r: _internalReactions)
        removeInternalReaction(r);
    
    
    
    ///Remove all species
    for(auto &m: _monomers) {
        for(auto &s : m->speciesFilamentVector()) _compartment->removeSpecies(s);
        for(auto &s : m->speciesBoundVector())    _compartment->removeSpecies(s);
        for(auto &s : m->speciesPlusEndVector())  _compartment->removeSpecies(s);
        for(auto &s : m->speciesMinusEndVector()) _compartment->removeSpecies(s);
    }
}


void CCylinder::addInternalReaction(ReactionBase* r) {
    //remove from compartment and chemsim
    _compartment->addInternalReactionUnique(std::unique_ptr<ReactionBase>(r));
    ChemSim::addReaction(ChemSimReactionKey(), r);
    
    ///add to local reaction list
    _internalReactions.push_back(r);
}

void CCylinder::addCrossCylinderReaction(CCylinder* other, ReactionBase* r) {
    //add to compartment and chemsim
    _compartment->addInternalReactionUnique(std::unique_ptr<ReactionBase>(r));
    ChemSim::addReaction(ChemSimReactionKey(), r);
    
    ///add to this reaction map
    _crossCylinderReactions[other].push_back(r);
    ///add to others reacting cylinders list
    other->addReactingCylinder(this);
}

void CCylinder::addReactingCylinder(CCylinder* other) { _reactingCylinders.push_back(other);}


void CCylinder::removeInternalReaction(ReactionBase* r) {
    ///remove from compartment and chemsim
    _compartment->removeInternalReaction(r);
    ChemSim::removeReaction(ChemSimReactionKey(), r);
    
    ///remove from internal reaction list
    auto it = std::find(_internalReactions.begin(), _internalReactions.end(), r);
    if (it != _internalReactions.end()) _internalReactions.erase(it);
}

void CCylinder:: removeAllInternalReactions() {
    
    for (auto &r : _internalReactions) {
        ///remove from compartment and chemsim
        _compartment->removeInternalReaction(r);
        ChemSim::removeReaction(ChemSimReactionKey(), r);
    }
    _internalReactions.clear();
}

void CCylinder::removeCrossCylinderReactions(CCylinder* other) {
    
    ///Remove from this map
    for(auto &r : _crossCylinderReactions[other]) {
        _compartment->removeInternalReaction(r);
        ChemSim::removeReaction(ChemSimReactionKey(), r);
    }
    _crossCylinderReactions.erase(other);
    
    ///also remove from reacting list of other ccylinder
    other->removeReactingCylinder(this);
}

void CCylinder::removeAllCrossCylinderReactions() {
    
    for(auto it = _crossCylinderReactions.begin(); it != _crossCylinderReactions.end(); it++) {
        
        ///Remove from this map
        for(auto &r : it->second) {
            _compartment->removeInternalReaction(r);
            ChemSim::removeReaction(ChemSimReactionKey(), r);
        }
        
        ///also remove from map of other ccylinder
        it->first->removeReactingCylinder(this);
    }
    _crossCylinderReactions.clear();
}

void CCylinder::removeReactingCylinder(CCylinder* other) {
    auto it = std::find(_reactingCylinders.begin(), _reactingCylinders.end(), other);
    if (it != _reactingCylinders.end()) _reactingCylinders.erase(it);
}


void CCylinder::updateReactions()
{
    ///loop through all reactions, passivate/activate
    for(auto &r : _internalReactions) {
        if(r->getProductOfReactants() == 0)
            r->passivateReaction();
        else
            r->activateReaction();
    }
    
    for(auto it = _crossCylinderReactions.begin(); it != _crossCylinderReactions.end(); it++) {
        
        for(auto &r : it->second) {
            if(r->getProductOfReactants() == 0)
                r->passivateReaction();
            else
                r->activateReaction();
        }
    }
}

void CCylinder::printCCylinder()
{
    std::cout << "Compartment:" << _compartment << std::endl;
    
    std::cout << "Composition of CCylinder: " << std::endl;
    for (auto &m : _monomers){
        m->print();
        std::cout << ":";
    }
}


