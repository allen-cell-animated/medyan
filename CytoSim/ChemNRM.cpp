//
//  ChemNRM.cpp
//  CytoSim
//
//  Created by Garegin Papoian on 5/10/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include <iostream>
#include "ChemNRM.h"
#include "ChemNRMImpl.h"
#include "utility.h"

ChemNRM::ChemNRM()
{
    _pimpl=make_unique<ChemNRMImpl>();
}

void ChemNRM::addReaction(Reaction *r){
    _pimpl->addReaction(r);
}

void ChemNRM::removeReaction(Reaction *r){
    _pimpl->removeReaction(r); 
}

void ChemNRM::run(int steps){
    _pimpl->run(steps);
}

size_t ChemNRM::getSize() const {
    return _pimpl->getSize();
}

float ChemNRM::getTime() const {
    return _pimpl->getTime();
}

void ChemNRM::initialize() {
    return _pimpl->initialize();
}

void ChemNRM::printReactions() const {
    return _pimpl->printReactions();
}
