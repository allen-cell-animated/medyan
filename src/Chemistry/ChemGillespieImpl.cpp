
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.0
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "ChemGillespieImpl.h"

RNodeGillespie::RNodeGillespie(ReactionBase *r, ChemGillespieImpl &chem_Gillespie)
    :_chem_Gillespie (chem_Gillespie), _react(r) {
    _react->setRnode(this);
    reset();
}

RNodeGillespie::~RNodeGillespie() noexcept {
    _react->setRnode(nullptr);
}


void RNodeGillespie::printSelf() const {
    cout << "RNodeGillespie: ptr=" << this << ", a=" << _a << ", a_penult=" <<
        _a_prev << ", points to Reaction:\n";
    cout << (*_react);
}

void RNodeGillespie::printDependents() const {
    cout << "RNodeGillespie: ptr=" << this <<
        ", the following RNodeGillespie objects are dependents:\n\n";
    for(auto rit = _react->dependents().begin();
        rit!=_react->dependents().end(); ++rit){
        
        RNodeGillespie *rn_other = (RNodeGillespie*)((*rit)->getRnode());
        rn_other->printSelf();
    }
    cout << endl;
}

void RNodeGillespie::activateReaction() {
    _chem_Gillespie.activateReaction(getReaction());
}

void RNodeGillespie::passivateReaction() {
    _chem_Gillespie.passivateReaction(getReaction());
}

void ChemGillespieImpl::initialize() {
    resetTime();
    for (auto &x : _map_rnodes){
        auto rn = x.second.get();
#ifdef TRACK_DEPENDENTS
        rn->getReaction()->activateReactionUnconditional();
#endif
        rn->reset();
    }
    _a_total = computeTotalA();
}


ChemGillespieImpl::~ChemGillespieImpl() noexcept{
    _map_rnodes.clear();
}

double ChemGillespieImpl::generateTau(double a) {
    exponential_distribution<double>::param_type pm(a);
    _exp_distr.param(pm);
    return _exp_distr(_eng);
}

double ChemGillespieImpl::generateUniform() {
    return _uniform_distr(_eng);
}

double ChemGillespieImpl::computeTotalA() {
    double rates_sum = 0;
    for (auto &x : _map_rnodes){
        auto rn = x.second.get();
        if(rn->getReaction()->isPassivated())
            continue;
        rates_sum+=rn->getPropensity();
    }
    return rates_sum;
}

bool ChemGillespieImpl::makeStep() {
    RNodeGillespie *rn_selected = nullptr;
    
    //Gillespie algorithm's first step; We assume that _a_total is up to date
    double tau = generateTau(_a_total);
    _t+=tau;
    syncGlobalTime();
    
    //Gillespie algorithm's second step: finding which reaction happened;
    double mu = _a_total*generateUniform();
    double rates_sum = 0;
    for (auto &x : _map_rnodes){
        
        auto rn = x.second.get();
        if(rn->getReaction()->isPassivated()) continue;
        
        rates_sum+=rn->getPropensity();
        
        if(rates_sum>mu){
            rn_selected = rn;
            break;
        }
    }
    if(rn_selected==nullptr){
        cout << "ChemGillespieImpl::makeStep() for loop: rates_sum=" <<
            rates_sum << ", mu=" << mu << ", _a_total=" << _a_total << endl;
        throw runtime_error("ChemGillespieImpl::makeStep(): No Reaction was selected during the Gillespie step!");
    }
    double a_new, a_penult;
    rn_selected->makeStep();
    if(!rn_selected->isPassivated()){
        rn_selected->reComputePropensity();
        a_new = rn_selected->getPropensity();
        a_penult = rn_selected->getPenultStepPropensity();
        _a_total = _a_total - a_penult + a_new;
    }

    // Updating dependencies
    ReactionBase *r = rn_selected->getReaction();
    
    if(r->updateDependencies()) {
        for(auto rit = r->dependents().begin(); rit!=r->dependents().end(); ++rit){
            RNodeGillespie *rn_other = (RNodeGillespie*)((*rit)->getRnode());
            rn_other->reComputePropensity();
            a_new = rn_other->getPropensity();
            a_penult = rn_other->getPenultStepPropensity();
            _a_total = _a_total - a_penult + a_new;
        }
    }
    
    // Send signal
#ifdef REACTION_SIGNALING
    r->emitSignal();
#endif
    syncGlobalTime();
    return true;
}

void ChemGillespieImpl::addReaction(ReactionBase *r) {
    _map_rnodes.emplace(r,make_unique<RNodeGillespie>(r,*this));
    ++_n_reacts;
}

void ChemGillespieImpl::removeReaction(ReactionBase *r) {
    _map_rnodes.erase(r);
    --_n_reacts;
}

void ChemGillespieImpl::printReactions() const {
    for (auto &x : _map_rnodes){
        auto rn = x.second.get();
        rn->printSelf();
    }
}

void ChemGillespieImpl::activateReaction(ReactionBase *r) {
    auto mit = _map_rnodes.find(r);
    if(mit!=_map_rnodes.end()){
        RNodeGillespie *rn_this = mit->second.get();
        rn_this->reComputePropensity();
        double a_new = rn_this->getPropensity();
        double a_penult = rn_this->getPenultStepPropensity();
        _a_total = _a_total - a_penult + a_new;
    }
    else
        throw out_of_range(
        "ChemGillespieImpl::activateReaction(...): Reaction not found!");
}

void ChemGillespieImpl::passivateReaction(ReactionBase *r) {
    auto mit = _map_rnodes.find(r);
    if(mit==_map_rnodes.end())
        throw out_of_range(
        "ChemGillespieImpl::passivateReaction(...): Reaction not found!");
    RNodeGillespie *rn_this = mit->second.get();
    
    double a_new, a_penult;
    a_penult = rn_this->getPropensity();
    rn_this->setPenultA(a_penult);
    a_new = 0;
    rn_this->setA(a_new);
    _a_total = _a_total - a_penult + a_new;
}
