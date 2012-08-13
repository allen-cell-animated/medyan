//
//  RSpecies.cpp
//  CytoSim
//
//  Created by Garegin Papoian on 5/22/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include <iostream>
#include "RSpecies.h"
#include "Reaction.h"

/// Print self into an iostream
std::ostream& operator<<(std::ostream& os, const chem::RSpecies& s){
    os << s.getFullName() << "[" << s.getN() << "]";
    return os;
}

namespace chem {
    
std::string RSpecies::getFullName() const {
    return _species.getFullName();
}


void RSpecies::activateAssocReactantReactions() {
    for (auto &r : _as_reactants)
        r->activateReaction();
}
    
void RSpecies::activateAssocProductReactions() {
    for (auto &r : _as_products)
        r->activateReaction();
}

void RSpecies::passivateAssocReactantReactions() {
    for (auto &r : _as_reactants)
        r->passivateReaction();
}
    
void RSpecies::passivateAssocProductReactions() {
    for (auto &r : _as_products)
        r->passivateReaction();
}
    
void RSpecies::startSignaling () {
    _signal = new RSpeciesCopyNChangedSignal;
}

void RSpecies::stopSignaling () {
    if (_signal!=nullptr)
        delete _signal;
    _signal = nullptr;
}
}
