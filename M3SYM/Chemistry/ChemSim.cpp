
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "ChemSim.h"

#include "ChemSimImpl.h"
    
ChemSimImpl* ChemSim::_pimpl = 0;
    
void ChemSim::setInstance(ChemSimImpl *csi) {
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