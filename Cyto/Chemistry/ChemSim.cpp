//
//  ChemSim.cpp
//  CytoSim
//
//  Created by Garegin Papoian on 5/10/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include "ChemSim.h"

#include "ChemSimImpl.h"
#include "utility.h"
    
ChemSimImpl* ChemSim::_pimpl = 0;
    
void ChemSim::setInstance(ChemSimImpl *csi)
{
    if (_pimpl != 0) delete _pimpl;
    _pimpl=csi;
}

void ChemSim::addReaction(ReactionBase *r){
    _pimpl->addReaction(r);
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

void ChemSim::printReactions() {
    return _pimpl->printReactions();
}