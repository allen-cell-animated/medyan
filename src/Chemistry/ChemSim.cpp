
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
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

bool ChemSim::run(double time){
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

//Dissipation
vector<double> ChemSim::getEnergy(){
    vector<double> returnvec;
    if(SysParams::Chemistry().dissTracking){
    returnvec.push_back(_pimpl->_dt->getCumDissEnergy());
    returnvec.push_back(_pimpl->_dt->getCumDissChemEnergy());
    returnvec.push_back(_pimpl->_dt->getCumDissMechEnergy());
    returnvec.push_back(_pimpl->_dt->getCumGChemEn());
    returnvec.push_back(_pimpl->_dt->getCumGMechEn());
    }
    return returnvec;
}


DissipationTracker* ChemSim::getDT(){
    return _pimpl->_dt;
}



