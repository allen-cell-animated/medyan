//
//  ChemSim.cpp
//  CytoSim
//
//  Created by Garegin Papoian on 5/10/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include <iostream>
#include "ChemSim.h"
#include "utility.h"

ChemSim::ChemSim(ChemSimImpl *csi)
{
    _pimpl=csi;
}

void ChemSim::addReaction(Reaction *r){
    _pimpl->addReaction(r);
}

void ChemSim::removeReaction(Reaction *r){
    _pimpl->removeReaction(r); 
}

void ChemSim::run(int steps){
    _pimpl->run(steps);
}

void ChemSim::initialize() {
    return _pimpl->initialize();
}

void ChemSim::printReactions() const {
    return _pimpl->printReactions();
}
