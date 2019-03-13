
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "ChemSim.h"

#include "ChemSimImpl.h"

void ChemSim::setInstance(ChemSimImpl *csi) {
    _pimpl=csi;
}

void ChemSim::addReaction(ReactionBase *r){
    _pimpl->addReaction(r);
}
    
void ChemSim::removeReaction(ReactionBase *r){
    _pimpl->removeReaction(r); 
}

bool ChemSim::run(floatingpoint time){
    return _pimpl->run(time);
}

bool ChemSim::runSteps(int steps){
    return _pimpl->runSteps(steps);
}

void ChemSim::initialize() {
    return _pimpl->initialize();
}

void ChemSim::printReactions() {
    return _pimpl->printReactions();
}
