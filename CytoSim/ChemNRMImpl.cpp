//
//  ChemNRMImpl.cpp
//  CytoSim
//
//  Created by Garegin Papoian on 5/6/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include <iostream>

//#define epsilon 

using namespace std;


#include "ChemNRMImpl.h"

RNodeNRM::RNodeNRM(Reaction *r, boost_heap &heap) : _react(r) {
    _react->setRnode(this);
    _handle = heap.emplace(this);
    _a = _react->computePropensity ();
}


void RNodeNRM::printSelf() const {
    cout << "RNodeNRM: ptr=" << this << ", tau=" << getTau() << ", a=" << _a << ", points to Reaction:\n";
    _react->printSelf();
}

void RNodeNRM::printDependents() const {
    cout << "RNodeNRM: ptr=" << this << ", the following RNodeNRM objects are dependents:\n\n";
    for(auto rit = _react->beginAffected(); rit!=_react->endAffected(); ++rit){
//        cout << "I am here [" << (*rit)->getRnode() << "]" << endl;
        RNodeNRM *rn_other = static_cast<RNodeNRM*>((*rit)->getRnode());
        rn_other->printSelf();
    }
    cout << endl;
}

void RNodeNRM::updateHeap(boost_heap &heap) {
    heap.update(_handle);
}

void ChemNRMImpl::initialize() {
    for (auto &x : _map_rnodes){
        auto rn = x.second.get();
        _generateNewRandTau(rn);
        rn->updateHeap(_heap);
    }
}

void ChemNRMImpl::_generateNewRandTau(RNodeNRM *rn) {
    rn->reComputePropensity();
    float a = rn->getPropensity();
    exponential_distribution<float>::param_type pm(a);
    _exp_distr.param(pm);
    float tau = _exp_distr(_eng) + _t;
    float prev_tau = rn->getTau();
    rn->setTau(tau);
//    cout << "ChemNRMImpl::generateNewRandTau(): RNodeNRM ptr=" << rn << ", was called. Prev tau=" << prev_tau 
//         << ", New tau=" << tau << endl;
}


void ChemNRMImpl::_makeStep() 
{
    RNodeNRM *rn = _heap.top()._rn;
    _t=rn->getTau();
    rn->makeStep();
    _generateNewRandTau(rn);
    rn->updateHeap(_heap);
//    cout << "ChemNRMImpl::makeStep(): RNodeNRM ptr=" << rn << " made a chemical step. " << endl;
//    rn->printSelf();
    // Now updating dependencies
    Reaction *r = rn->getReaction();
    for(auto rit = r->beginAffected(); rit!=r->endAffected(); ++rit){
        RNodeNRM *rn_other = static_cast<RNodeNRM*>((*rit)->getRnode());
        float a_old = rn_other->getPropensity();
        rn_other->reComputePropensity();
        float tau_old = rn_other->getTau();
        float a_new = rn_other->getPropensity();
        float tau_new = (a_old/a_new)*(tau_old-_t)+_t; 
        rn_other->setTau(tau_new);
        rn_other->updateHeap(_heap);
    }
}

void ChemNRMImpl::printReactions() const {
    for (auto &x : _map_rnodes){
        auto rn = x.second.get();
        rn->printSelf();
    }
}

