//
//  ChemSim.cpp
//  CytoSim
//
//  Created by Garegin Papoian on 5/10/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include <iostream>
#include "ChemSim.h"
#include "ChemSimImpl.h"
#include "utility.h"

namespace chem {
    
ChemSim::ChemSim(ChemSimImpl *csi)
{
    _pimpl=csi;
}

void ChemSim::addReaction(ReactionBase *r){
    _pimpl->addReaction(r);
}
    
void ChemSim::addAndActivateReaction(ReactionBase* r){
    _pimpl->addAndActivateReaction(r);
}
    
void ChemSim::removeReaction(ReactionBase *r){
    _pimpl->removeReaction(r); 
}

bool ChemSim::run(int steps){
    return _pimpl->run(steps);
}

void ChemSim::initialize() {
    return _pimpl->initialize();
}

void ChemSim::printReactions() const {
    return _pimpl->printReactions();
}

} // end of namespace