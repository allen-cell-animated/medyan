//
//  ChemNRMImpl.cpp
//  CytoSim
//
//  Created by Garegin Papoian on 5/6/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include <iostream>
#include <chrono>

//#define epsilon 

using namespace std;


#include "ChemNRMImpl.h"

namespace chem {
    
    long long int elapsed14 = 0;
    long long int elapsed16 = 0;
    long long int elapsed23 = 0;
    long long int elapsed45 = 0;
    long long int elapsed56 = 0;

RNodeNRM::RNodeNRM(Reaction *r, ChemNRMImpl &chem_nrm) :
 _chem_nrm (chem_nrm), _react(r) {
    _react->setRnode(this);
     boost_heap *heap = _chem_nrm.getHeap();
    _handle = heap->emplace(this);
    _a = _react->computePropensity ();
}

RNodeNRM::~RNodeNRM () 
{
//    cout << "RNodeNRM::~RNodeNRM (): ptr=" << this << ", tau=" << getTau() << ", a=" << _a << ", points to Reaction:\n";
//    _react->printSelf();
    boost_heap *heap = _chem_nrm.getHeap();
    heap->erase(_handle);
    _react->setRnode(nullptr);
}


void RNodeNRM::printSelf() const {
    cout << "RNodeNRM: ptr=" << this << ", tau=" << getTau() << ", a=" << _a << ", points to Reaction:\n";
    cout << (*_react);
}

void RNodeNRM::printDependents() const {
    cout << "RNodeNRM: ptr=" << this << ", the following RNodeNRM objects are dependents:\n\n";
    for(auto rit = _react->dependents().begin(); rit!=_react->dependents().end(); ++rit){
//        cout << "I am here [" << (*rit)->getRnode() << "]" << endl;
        RNodeNRM *rn_other = static_cast<RNodeNRM*>((*rit)->getRnode());
        rn_other->printSelf();
    }
    cout << endl;
}

void RNodeNRM::updateHeap() {
    boost_heap *heap = _chem_nrm.getHeap();
    heap->update(_handle);
}

void RNodeNRM::generateNewRandTau() {
    double tau;
    reComputePropensity();//calculated new _a
#ifdef TRACK_ZERO_COPY_N
    tau = _chem_nrm.generateTau(_a) + _chem_nrm.getTime();
#else
    if(_a<1.0e-10) // std::numeric_limits< double >::min()
        tau = numeric_limits<double>::infinity();
    else
        tau = _chem_nrm.generateTau(_a) + _chem_nrm.getTime();
#endif
    setTau(tau);
    //    cout << "RNodeNRM::generateNewRandTau(): RNodeNRM ptr=" << this << ", was called. Prev tau=" << prev_tau 
    //         << ", New tau=" << tau << endl;
}

void RNodeNRM::activateReaction() {
    generateNewRandTau();
    updateHeap();
//    cout << "RNodeNRM::activateReaction(): ptr=" << _react << ", has been activated." << endl; 
//    cout << (*_react);
//    cout << "EndOfRNodeNRM::activateReaction..." << endl;

}

void RNodeNRM::passivateReaction() {
    _a=0;
    double tau = numeric_limits<double>::infinity();
    setTau(tau);
    updateHeap();
//    cout << "RNodeNRM::passivateReaction(): ptr=" << _react << ", has been passivated, and tau set to inf." << endl; 
//    cout << (*_react);
//    cout << "EndOfRNodeNRM::passivateReaction..." << endl;
}

void ChemNRMImpl::initialize() {
    resetTime();
    for (auto &x : _map_rnodes){
        auto rn = x.second.get();
        rn->getReaction()->activateReaction();
    }
}


ChemNRMImpl::~ChemNRMImpl() {
    _map_rnodes.clear();
//    double sec = std::pow(10,9);
//    cout << "Elapsed times in ChemNRMImpl::makeStep(): dt14=" << elapsed14/sec << ", dt23=" << elapsed23/sec  << ", dt45=" << elapsed45/sec <<  ", dt56=" << elapsed56/sec <<  ", dt16=" << elapsed16/sec << endl;
}

double ChemNRMImpl::generateTau(double a){
    exponential_distribution<double>::param_type pm(a);
    _exp_distr.param(pm);
    return _exp_distr(_eng);
}

bool ChemNRMImpl::makeStep()
{
//    cout << "[ChemNRMImpl::_makeStep(): Starting..." << endl;
//    std::chrono::time_point<std::chrono::system_clock> chk1, chk2, chk3, chk4, chk5, chk6;

//    chk1 = std::chrono::system_clock::now();

    RNodeNRM *rn = _heap.top()._rn;
    double tau_top = rn->getTau();
    if(tau_top==numeric_limits<double>::infinity()){
        cout << "The heap has been exhausted - no more reactions to fire, returning..." << endl;
        return false;
    }
    _t=tau_top;
    rn->makeStep();
    if(!rn->isPassivated()){
        rn->generateNewRandTau();
        rn->updateHeap();
    }
    
//    chk4 = std::chrono::system_clock::now();
    
//    cout << "ChemNRMImpl::makeStep(): RNodeNRM ptr=" << rn << " made a chemical step. t=" << _t << "\n" << endl;
//    rn->printSelf();
    // Updating dependencies
    Reaction *r = rn->getReaction();
    for(auto rit = r->dependents().begin(); rit!=r->dependents().end(); ++rit){
        RNodeNRM *rn_other = static_cast<RNodeNRM*>((*rit)->getRnode());
        double a_old = rn_other->getPropensity();
//        chk2 = std::chrono::system_clock::now();
        rn_other->reComputePropensity();
//        chk3 = std::chrono::system_clock::now();
        double tau_new;
        double tau_old = rn_other->getTau();
        double a_new = rn_other->getPropensity();
#ifdef TRACK_ZERO_COPY_N
        tau_new = (a_old/a_new)*(tau_old-_t)+_t;
#else
        if(a_new<1.0e-15) // std::numeric_limits< double >::min()
            tau_new = numeric_limits<double>::infinity();
        else if (a_old<1.0e-15){
            rn_other->generateNewRandTau();
            tau_new = rn_other->getTau();
        }
        else{
            tau_new = (a_old/a_new)*(tau_old-_t)+_t;
        }
#endif
        rn_other->setTau(tau_new);
        rn_other->updateHeap();
//        elapsed23 += std::chrono::duration_cast<std::chrono::nanoseconds>(chk3-chk2).count();
    }
    
//    chk5 = std::chrono::system_clock::now();
    
    // Send signals
    r->emitSignal();
#ifdef RSPECIES_SIGNALING
    for(auto sit = r->beginReactants(); sit!=r->endReactants(); ++sit){
        if((*sit)->isSignaling())
            (*sit)->emitSignal(-1);
    }
    for(auto sit = r->beginProducts(); sit!=r->endProducts(); ++sit){
        if((*sit)->isSignaling())
            (*sit)->emitSignal(1);
    }
#endif

//    cout << "ChemNRMImpl::_makeStep(): Ending...]\n\n" << endl;
    syncGlobalTime();
    
//    chk6 = std::chrono::system_clock::now();
    
//    elapsed14 += std::chrono::duration_cast<std::chrono::nanoseconds>(chk4-chk1).count();
//    elapsed45 += std::chrono::duration_cast<std::chrono::nanoseconds>(chk5-chk4).count();
//    elapsed56 += std::chrono::duration_cast<std::chrono::nanoseconds>(chk6-chk5).count();
//    elapsed16 += std::chrono::duration_cast<std::chrono::nanoseconds>(chk6-chk1).count();

    return true;
}

void ChemNRMImpl::addReaction(Reaction *r) {
    _map_rnodes.emplace(r,make_unique<RNodeNRM>(r,*this));
    ++_n_reacts;
}

void ChemNRMImpl::removeReaction(Reaction *r) {
    _map_rnodes.erase(r);
    --_n_reacts;
}

void ChemNRMImpl::printReactions() const {
    for (auto &x : _map_rnodes){
        auto rn = x.second.get();
        rn->printSelf();
    }
}

} // end of namespace