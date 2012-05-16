//
//  Species.cpp
//  CytoSim-Experimenting
//
//  Created by Garegin Papoian on 4/21/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include <iostream>

#include "Species.h"
#include "Reaction.h"
#include "Signaling.h"

std::vector<std::string> SpeciesType::_vec_type_name = {"Bulk", "Diffusing", "Membrane", "Filament", "Walking", "Motors" };

void Species::activateAssocReactions() {
    for (auto &r : _as_reactants)
        r->activateReaction();
}

void Species::passivateAssocReacts() {
    for (auto &r : _as_reactants)
        r->passivateReaction();
}

void Species::makeSignaling (SignalingManager &sm) {
    sm.addSignalingSpecies(this);
    _is_signaling=true;
}

void Species::stopSignaling (SignalingManager &sm) {
    sm.disconnect_semiprivate(this);
    _is_signaling=false;
}

void Species::printSelf() const {
    std::cout << "Species: " << getFullName() << ", copy_num=" << (int)this->_n << ", ptr=" << (void*)this << std::endl;
//    std::cout << "This species is involved in the following Reactions: \nAs a reactant:\n";
//    for(auto &r : _freactions)
//        std::cout << r << ",";
//    std::cout << "\nAs a product:\n";
//    for(auto &r : _breactions)
//        std::cout << r << ",";
//    std::cout << "\n";
}

