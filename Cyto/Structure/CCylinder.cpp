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
    for(auto &r: rhs._reactions)
        addReaction(r->clone(c->speciesContainer()));
    for(auto &r: rhs._frontReactions)
        addFrontReaction(r->clone(c->speciesContainer()), true);
    for(auto &r: rhs._backReactions)
        addBackReaction(r->clone(c->speciesContainer()), true);
    
    ///Update and return
    this->updateReactions();
}

CCylinder::~CCylinder()
{
    ///Remove all reactions
    for(auto &r: _reactions)
        removeReaction(r);
    for(auto &r: _frontReactions)
        removeReaction(r);
    for(auto &r: _backReactions)
        removeReaction(r);
    
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


void CCylinder::addReaction(ReactionBase* r) {
    //remove from compartment and chemsim
    _compartment->addInternalReactionUnique(std::unique_ptr<ReactionBase>(r));
    ChemSim::addReaction(ChemSimReactionKey(), r);
    
    ///add to local reaction list
    _reactions.push_back(r);
}

void CCylinder::addFrontReaction(ReactionBase* r, bool manage) {
    //remove from compartment and chemsim
    if(manage) {
        _compartment->addInternalReactionUnique(std::unique_ptr<ReactionBase>(r));
        ChemSim::addReaction(ChemSimReactionKey(), r);
    }
    
    ///add to local reaction list
    _frontReactions.push_back(r);
}

void CCylinder::addBackReaction(ReactionBase* r, bool manage) {
    //remove from compartment and chemsim
    if(manage) {
        _compartment->addInternalReactionUnique(std::unique_ptr<ReactionBase>(r));
        ChemSim::addReaction(ChemSimReactionKey(), r);
    }
    
    ///add to local reaction list
    _backReactions.push_back(r);
}

void CCylinder::removeReaction(ReactionBase* r) {
    ///remove from compartment and chemsim
    _compartment->removeInternalReaction(r);
    ChemSim::removeReaction(ChemSimReactionKey(), r);
}

void CCylinder::updateReactions()
{
    ///loop through all reactions, passivate/activate
    for(auto &r : _backReactions) {
        if(r->getProductOfReactants() == 0)
            r->passivateReaction();
        else
            r->activateReaction();
    }
    for(auto &r : _reactions) {
        if(r->getProductOfReactants() == 0)
            r->passivateReaction();
        else
            r->activateReaction();
    }
    for(auto &r : _frontReactions) {
        if(r->getProductOfReactants() == 0)
            r->passivateReaction();
        else
            r->activateReaction();
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


