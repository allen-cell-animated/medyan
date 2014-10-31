//
//  CCylinder.cpp
//  Cyto
//
//  Created by James Komianos on 9/22/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "CCylinder.h"

CCylinder::CCylinder(const CCylinder& rhs, Compartment* c) : _compartment(c), _pCylinder(rhs._pCylinder)
{
    ///copy all monomers, bounds
    for(auto &m : rhs._monomers)
        _monomers.push_back(std::unique_ptr<CMonomer>(m->clone(c)));
    
    ///copy all reactions
    for(auto &r: rhs._internalReactions) addInternalReaction(r->clone(c->speciesContainer()));
    for(auto it = rhs._crossCylinderReactions.begin(); it != rhs._crossCylinderReactions.end(); it++) {
        ///Copy map
        for(auto &r : it->second) addCrossCylinderReaction(it->first, r->clone(c->speciesContainer()));
    }
    
    ///Copy reacting cylinders, Clone reactions where this cylinder is involved
    for(auto &ccyl : rhs._reactingCylinders) {
        ///add as a reacting cylinder
        addReactingCylinder(ccyl);
        ///clone reactions
        for(auto &r: ccyl->getCrossCylinderReactions()[const_cast<CCylinder*>(&rhs)])
            ccyl->addCrossCylinderReaction(this, r->clone(c->speciesContainer()));
    }
    
    ///Update and return
    this->updateReactions();
}

CCylinder::~CCylinder()
{
    ///remove these reactions as dependents to other ccylinders
    passivateReactions();
    
    ///Remove all reactions owned by this ccylinder
    removeAllInternalReactions();
    removeAllCrossCylinderReactions();
    
    ///remove all reactions involving this ccylinder
    removeAllReactingCylinders();
    
    ///Remove all species
    for(auto &m: _monomers) {
        
        short numFilamentSpecies = SystemParameters::Chemistry().numFilamentSpecies;
        for(int i = 0; i < numFilamentSpecies; i++) {
            SpeciesFilament* s = m->speciesFilament(i);
            if(s != nullptr) _compartment->removeSpecies(s);
        }
        short numPlusEndSpecies = SystemParameters::Chemistry().numPlusEndSpecies;
        for(int i = 0; i < numPlusEndSpecies; i++) {
            SpeciesPlusEnd* s = m->speciesPlusEnd(i);
            if(s != nullptr) _compartment->removeSpecies(s);
        }
        short numMinusEndSpecies = SystemParameters::Chemistry().numMinusEndSpecies;
        for(int i = 0; i < numMinusEndSpecies; i++) {
            SpeciesMinusEnd* s = m->speciesMinusEnd(i);
            if(s != nullptr) _compartment->removeSpecies(s);
        }
        short numBoundSpecies = SystemParameters::Chemistry().numBoundSpecies;
        for(int i = 0; i < numBoundSpecies; i++) {
            SpeciesBound* s = m->speciesBound(i);
            if(s != nullptr) _compartment->removeSpecies(s);
        }
        short numLinkerSpecies = SystemParameters::Chemistry().numLinkerSpecies;
        for(int i = 0; i < numLinkerSpecies; i++) {
            SpeciesLinker* s = m->speciesLinker(i);
            if(s != nullptr) _compartment->removeSpecies(s);
        }
        short numMotorSpecies = SystemParameters::Chemistry().numMotorSpecies;
        for(int i = 0; i < numMotorSpecies; i++) {
            SpeciesMotor* s = m->speciesMotor(i);
            if(s != nullptr) _compartment->removeSpecies(s);
        }
    }
}


void CCylinder::addInternalReaction(ReactionBase* r) {
    //remove from compartment and chemsim
    _compartment->addInternalReactionUnique(std::unique_ptr<ReactionBase>(r));
    ChemSim::addReaction(ChemSimReactionKey(), r);
    
    ///add to local reaction list
    _internalReactions.insert(r);
}

void CCylinder::addCrossCylinderReaction(CCylinder* other, ReactionBase* r) {
    //add to compartment and chemsim
    _compartment->addInternalReactionUnique(std::unique_ptr<ReactionBase>(r));
    ChemSim::addReaction(ChemSimReactionKey(), r);
    
    ///add to this reaction map
    _crossCylinderReactions[other].insert(r);
    ///add to others reacting cylinders list
    other->addReactingCylinder(this);
}

void CCylinder::addReactingCylinder(CCylinder* other) { _reactingCylinders.insert(other);}


void CCylinder::removeInternalReaction(ReactionBase* r) {
    ///remove from compartment and chemsim
    _compartment->removeInternalReaction(r);
    ChemSim::removeReaction(ChemSimReactionKey(), r);
    
    ///remove from internal reaction list
    auto it = _internalReactions.find(r);
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

void CCylinder::removeCrossCylinderReaction(CCylinder* other, ReactionBase* r) {

    auto it = _crossCylinderReactions[other].find(r);
    if(it != _crossCylinderReactions[other].end()) _crossCylinderReactions[other].erase(it);
    
    _compartment->removeInternalReaction(r);
    ChemSim::removeReaction(ChemSimReactionKey(), r);
    
}


void CCylinder::removeCrossCylinderReactions(CCylinder* other) {
    
    ///Remove from this map
    for(auto &r : _crossCylinderReactions[other]) {
        _compartment->removeInternalReaction(r);
        ChemSim::removeReaction(ChemSimReactionKey(), r);
    }
    _crossCylinderReactions.erase(other);
    
    ///also remove from reacting list of other ccylinder
    auto it = other->_reactingCylinders.find(this);
    if(it != other->_reactingCylinders.end()) other->_reactingCylinders.erase(it);
}

void CCylinder::removeAllCrossCylinderReactions() {
    
    for(auto it = _crossCylinderReactions.begin(); it != _crossCylinderReactions.end(); it++) {
        
        ///Remove from this map
        for(auto &r : it->second) {
            _compartment->removeInternalReaction(r);
            ChemSim::removeReaction(ChemSimReactionKey(), r);
        }
        
        ///also remove from list of other ccylinder
        auto it2 = it->first->_reactingCylinders.find(this);
        if(it2 != it->first->_reactingCylinders.end()) it->first->_reactingCylinders.erase(it2);
    }
    _crossCylinderReactions.clear();
}

void CCylinder::removeReactingCylinder(CCylinder* other) { other->removeCrossCylinderReactions(this);}

void CCylinder::removeAllReactingCylinders() {
    auto tempReactingCylinders = _reactingCylinders;
    for(auto &cc : tempReactingCylinders) cc->removeCrossCylinderReactions(this);
}


void CCylinder::updateReactions()
{
    ///loop through all reactions, passivate/activate
    for(auto &r : _internalReactions) {
        if(r->getProductOfReactants() == 0) r->passivateReaction();
        else r->activateReaction();
    }
    
    for(auto it = _crossCylinderReactions.begin(); it != _crossCylinderReactions.end(); it++) {
        
        for(auto &r : it->second) {
            if(r->getProductOfReactants() == 0) r->passivateReaction();
            else r->activateReaction();
        }
    }
    
    for(auto &cc : _reactingCylinders) {
        
        for (auto &r : cc->_crossCylinderReactions[this]) {
            if(r->getProductOfReactants() == 0) r->passivateReaction();
            else r->activateReaction();
        }
    }
}

void CCylinder::passivateReactions()
{
    ///loop through all reactions, passivate/activate
    for(auto &r : _internalReactions) r->passivateReaction();
    
    for(auto it = _crossCylinderReactions.begin(); it != _crossCylinderReactions.end(); it++)
        for(auto &r : it->second) r->passivateReaction();
    
    for(auto &cc : _reactingCylinders)
        for (auto &r : cc->_crossCylinderReactions[this]) r->passivateReaction();
    
}

void CCylinder::printCCylinder()
{
    std::cout << "Compartment:" << _compartment << std::endl;
    
    std::cout << "Composition of CCylinder: " << std::endl;
    for (auto &m : _monomers){
        m->print();
        std::cout << ":";
    }
    std::cout << std::endl;
}


