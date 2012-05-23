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
#include "Signaling.h"

/// Print self into an iostream
std::ostream& operator<<(std::ostream& os, const chem::RSpecies& s){
    os << s.getFullName() << "[" << s.getN() << "]";
    return os;
}

namespace chem {
    
std::string RSpecies::getFullName() const {
    return _species.getFullName();
}


void RSpecies::activateAssocReactions() {
    for (auto &r : _as_reactants)
        r->activateReaction();
}

void RSpecies::passivateAssocReacts() {
    for (auto &r : _as_reactants)
        r->passivateReaction();
}

void RSpecies::makeSignaling (ChemSignal &sm) {
    sm.addSignalingRSpecies(this);
    _is_signaling=true;
}

void RSpecies::stopSignaling (ChemSignal &sm) {
    sm.disconnect(this);
    _is_signaling=false;
}
}
