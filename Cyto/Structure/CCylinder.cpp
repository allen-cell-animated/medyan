//
//  CCylinder.cpp
//  Cyto
//
//  Created by James Komianos on 9/22/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "CCylinder.h"

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
    
    ///remove from local list
    _reactions.erase(std::find(_reactions.begin(), _reactions.end(), r));
}

void CCylinder::removeFrontReaction(ReactionBase* r, bool manage) {
    ///remove from compartment and chemsim
    if(manage) {
        _compartment->removeInternalReaction(r);
        ChemSim::removeReaction(ChemSimReactionKey(), r);
    }
    
    ///remove from local list
    _frontReactions.erase(std::find(_reactions.begin(), _reactions.end(), r));
}

void CCylinder::removeBackReaction(ReactionBase* r, bool manage) {
    ///remove from compartment and chemsim
    if(manage) {
        _compartment->removeInternalReaction(r);
        ChemSim::removeReaction(ChemSimReactionKey(), r);
    }
    
    ///remove from local list
    _backReactions.erase(std::find(_reactions.begin(), _reactions.end(), r));
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


