
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

#include <iostream>
#include <chrono>
#include "Rand.h"
#ifdef BOOST_MEM_POOL
    #include <boost/pool/pool.hpp>
    #include <boost/pool/pool_alloc.hpp>
    #include <boost/math/special_functions/fpclassify.hpp>
#endif

#include "ChemNRMImpl.h"

#ifdef BOOST_MEM_POOL
#ifdef BOOST_POOL_MEM_RNODENRM
boost::pool<> allocator_rnodenrm(sizeof(RNodeNRM),BOOL_POOL_NSIZE);
#endif
    
#ifdef BOOST_POOL_MEM_PQNODE
boost::pool<> allocator_pqnode(sizeof(PQNode),BOOL_POOL_NSIZE);
#endif


#ifdef BOOST_POOL_MEM_PQNODE
void* PQNode::operator new(size_t size) {
    void *ptr = boost::fast_pool_allocator<PQNode>::allocate();
    return ptr;
}

void PQNode::operator delete(void* ptr) noexcept {
    boost::fast_pool_allocator<PQNode>::deallocate((PQNode*)ptr);
}
#endif
 
#ifdef BOOST_POOL_MEM_RNODENRM
void* RNodeNRM::operator new(size_t size) {
    void *ptr = boost::fast_pool_allocator<RNodeNRM>::allocate();
    return ptr;
}

void RNodeNRM::operator delete(void* ptr) noexcept {
    boost::fast_pool_allocator<RNodeNRM>::deallocate((RNodeNRM*)ptr);
}
#endif
#endif

RNodeNRM::RNodeNRM(ReactionBase *r, ChemNRMImpl &chem_nrm)
    : _chem_nrm (chem_nrm), _react(r) {
    _react->setRnode(this);
     boost_heap *heap = _chem_nrm.getHeap();
    _handle = heap->emplace(this);
    _a = _react->computePropensity();
}

RNodeNRM::~RNodeNRM() noexcept {
    boost_heap *heap = _chem_nrm.getHeap();
    heap->erase(_handle);
    _react->setRnode(nullptr);
}

void RNodeNRM::printSelf() const {
    cout << "RNodeNRM: ptr=" << this << ", tau=" << getTau() <<
        ", a=" << _a << ", points to Reaction:\n";
    cout << (*_react);
}

void RNodeNRM::printDependents() const {
    cout << "RNodeNRM: ptr=" << this
    << ", the following RNodeNRM objects are dependents:\n\n";
    for(auto rit = _react->dependents().begin();
        rit!=_react->dependents().end(); ++rit){
        RNodeNRM *rn_other = (RNodeNRM*)((*rit)->getRnode());
        rn_other->printSelf();
    }
    cout << endl;
}

void RNodeNRM::updateHeap() {
    boost_heap *heap = _chem_nrm.getHeap();
    heap->update(_handle);
}

void RNodeNRM::generateNewRandTau() {
    double newTau;
    reComputePropensity();//calculated new _a
    
#ifdef TRACK_ZERO_COPY_N
    newTau = _chem_nrm.generateTau(_a) + _chem_nrm.getTime();
#else
    if(_a<1.0e-10) // numeric_limits< double >::min()
        newTau = numeric_limits<double>::infinity();
    else
        newTau = _chem_nrm.generateTau(_a) + _chem_nrm.getTime();
#endif
    setTau(newTau);
}

void RNodeNRM::activateReaction() {
    generateNewRandTau();
    updateHeap();
}

void RNodeNRM::passivateReaction() {
    _a=0;
    double tau = numeric_limits<double>::infinity();
    setTau(tau);
    updateHeap();
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
}

double ChemNRMImpl::generateTau(double a){
    exponential_distribution<double>::param_type pm(a);
    
    _exp_distr.param(pm);
    return _exp_distr(Rand::_eng);
}

bool ChemNRMImpl::makeStep() {

    //try to get a reaction
    if(_heap.empty()) {
        cout << "There are no reactions to fire, returning..." << endl;
        return false;
    }
    RNodeNRM *rn = _heap.top()._rn;
    double tau_top = rn->getTau();
    if(tau_top==numeric_limits<double>::infinity()){
        cout << "The heap has been exhausted - no more reactions to fire, returning..." << endl;
        return false;
    }    
    ///DEBUG
    //assert heap ordering
    if(tau_top < _t) {
        cout << "ERROR: The heap is not correctly sorted, returning..." << endl;
        cout << "Tau top = " << tau_top << endl;
        cout << "Tau current = " << _t << endl;
        cout << "Reaction type = " << rn->getReaction()->getReactionType() << endl;
        rn->printSelf();
        return false;
    }
    
    double t_prev = _t;
    
    _t=tau_top;
    syncGlobalTime();
    
    rn->makeStep();
#if defined TRACK_ZERO_COPY_N || defined TRACK_UPPER_COPY_N
    if(!rn->isPassivated()){
#endif
        rn->generateNewRandTau();
        rn->updateHeap();
#if defined TRACK_ZERO_COPY_N || defined TRACK_UPPER_COPY_N
    }
#endif
    // Updating dependencies
    ReactionBase *r = rn->getReaction();
    
    if(r->updateDependencies()) {
    
        for(auto rit = r->dependents().begin(); rit!=r->dependents().end(); ++rit){
            
            RNodeNRM *rn_other = (RNodeNRM*)((*rit)->getRnode());
            double a_old = rn_other->getPropensity();
            
            //recompute propensity
            rn_other->reComputePropensity();
            
            double tau_new;
            double tau_old = rn_other->getTau();
            
            double a_new = rn_other->getPropensity();
            
#ifdef TRACK_ZERO_COPY_N
            //recompute tau
            tau_new = (a_old/a_new)*(tau_old-_t)+_t;
#else
            //recompute tau
            if(a_new<1.0e-15) // numeric_limits< double >::min()
                tau_new = numeric_limits<double>::infinity();
            else if (a_old<1.0e-15){
                rn_other->generateNewRandTau();
                tau_new = rn_other->getTau();
            }
            else{
                tau_new = (a_old/a_new)*(tau_old-_t)+_t;
            }
#endif
            if(boost::math::isnan(tau_new)){tau_new=numeric_limits<double>::infinity();}
            ///DEBUG
            if(tau_new < _t) {
                
                cout << "WARNING: Generated tau may be incorrect. " << endl;
                
                cout << "Tau new = " << tau_new << endl;
                cout << "Tau old = " << tau_old << endl;
                cout << "Current global t = " << _t << endl;
                cout << "Previous global t = " << t_prev << endl;
                cout << "a_old = " << a_old << endl;
                cout << "a_new = " << a_new << endl;
                
                cout << "Reaction type = " << rn->getReaction()->getReactionType() << endl;
                
                
                rn->printSelf();
                rn_other->printSelf();
            }
            
            rn_other->setTau(tau_new);
            rn_other->updateHeap();
        }
    }
#ifdef REACTION_SIGNALING
    // Send signal
    r->emitSignal();
#endif
    return true;
}

void ChemNRMImpl::addReaction(ReactionBase *r) {
    _map_rnodes.emplace(r,make_unique<RNodeNRM>(r,*this));
    ++_n_reacts;
}

void ChemNRMImpl::removeReaction(ReactionBase *r) {
    _map_rnodes.erase(r);
    --_n_reacts;
}

void ChemNRMImpl::printReactions() const {
    for (auto &x : _map_rnodes){
        auto rn = x.second.get();
        rn->printSelf();
    }
}
