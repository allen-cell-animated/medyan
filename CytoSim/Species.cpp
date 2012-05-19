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

namespace chem {
    
void Species::activateAssocReactions() {
    for (auto &r : _as_reactants)
        r->activateReaction();
}

void Species::passivateAssocReacts() {
    for (auto &r : _as_reactants)
        r->passivateReaction();
}

void Species::makeSignaling (ChemSignal &sm) {
    sm.addSignalingSpecies(this);
    _is_signaling=true;
}

void Species::stopSignaling (ChemSignal &sm) {
    sm.disconnect_semiprivate(this);
    _is_signaling=false;
}

}
