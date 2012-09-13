//
//  ChemGillespieImpl.cpp
//  CytoSim
//
//  Created by Garegin Papoian on 8/11/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include "ChemGillespieImpl.h"

#include <iostream>

using namespace std;

namespace chem {
    
    RNodeGillespie::RNodeGillespie(ReactionBase *r, ChemGillespieImpl &chem_Gillespie) :
    _chem_Gillespie (chem_Gillespie), _react(r) {
        _react->setRnode(this);
        reset();
    }
    
    RNodeGillespie::~RNodeGillespie ()
    {
        //    cout << "RNodeGillespie::~RNodeGillespie (): ptr=" << this << ", tau=" << getTau() << ", a=" << _a << ", points to Reaction:\n";
        //    _react->printSelf();
        _react->setRnode(nullptr);
    }
    
    
    void RNodeGillespie::printSelf() const {
        cout << "RNodeGillespie: ptr=" << this << ", a=" << _a << ", a_penult=" << _a_prev << ", points to Reaction:\n";
        cout << (*_react);
    }
    
    void RNodeGillespie::printDependents() const {
        cout << "RNodeGillespie: ptr=" << this << ", the following RNodeGillespie objects are dependents:\n\n";
        for(auto rit = _react->dependents().begin(); rit!=_react->dependents().end(); ++rit){
            //        cout << "RNodeGillespie::printDependents():" << (*rit)->getRnode() << "]" << endl;
            RNodeGillespie *rn_other = static_cast<RNodeGillespie*>((*rit)->getRnode());
            rn_other->printSelf();
        }
        cout << endl;
    }
    
    void RNodeGillespie::activateReaction() {
        _chem_Gillespie.activateReaction(getReaction());
        // reComputePropensity();
        //    cout << "RNodeGillespie::activateReaction(): ptr=" << _react << ", has been activated." << endl;
        //    cout << (*_react);
        //    cout << "EndOfRNodeGillespie::activateReaction..." << endl;
        
    }
    
    void RNodeGillespie::passivateReaction() {
        _chem_Gillespie.passivateReaction(getReaction());
        //_a=0;
        //    cout << "RNodeGillespie::passivateReaction(): ptr=" << _react << ", has been passivated, and tau set to inf." << endl;
        //    cout << (*_react);
        //    cout << "EndOfRNodeGillespie::passivateReaction..." << endl;
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
    
    
    ChemGillespieImpl::~ChemGillespieImpl() {
        _map_rnodes.clear();
    }
    
    double ChemGillespieImpl::generateTau(double a) {
        exponential_distribution<double>::param_type pm(a);
        _exp_distr.param(pm);
        return _exp_distr(_eng);
    }
    
    double ChemGillespieImpl::generateUniform()  {
        return _uniform_distr(_eng);
    }
    
    double ChemGillespieImpl::computeTotalA(){
        double rates_sum = 0;
        for (auto &x : _map_rnodes){
            auto rn = x.second.get();
            if(rn->getReaction()->isPassivated())
                continue;
            rates_sum+=rn->getPropensity();
        }
        return rates_sum;
    }
    
    bool ChemGillespieImpl::makeStep()
    {
//        cout << "\n\n[ChemGillespieImpl::_makeStep(): Starting..." << ", _a_total=" << _a_total << endl;
//        printReactions();
        
//        if(_a_total!=computeTotalA()){
//            cout << "ChemGillespieImpl::makeStep(): " << _a_total << " vs " << computeTotalA() << endl;
//            assert(0 && "ChemGillespieImpl::makeStep(): The current total rate is not consistent");
//        }
        
        //_a_total=computeTotalA();
        
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
            if(rn->getReaction()->isPassivated())
                continue;
            rates_sum+=rn->getPropensity();
//            rn->printSelf();
//            cout << rn->getReaction()->computePropensity() << " " << rn->getPropensity() << endl;
//            if(abs(rn->getReaction()->computePropensity()-rn->getPropensity())>0.00001){
//                rn->printSelf();
//                cout << "Problem: " << rn->getReaction()->computePropensity() << " " << rn->getPropensity() << endl;
//            }
//            cout << "ChemGillespieImpl::makeStep() for loop: a=" << rn->getPropensity() << ", rates_sum=" << rates_sum << ", mu="
//                << mu << ", _a_total=" << _a_total << endl;
            if(rates_sum>mu){
                rn_selected = rn;
                break;
            }
        }
                
        if(rn_selected==nullptr){
            cout << "ChemGillespieImpl::makeStep() for loop: rates_sum=" << rates_sum << ", mu="
                        << mu << ", _a_total=" << _a_total << endl;
            throw std::runtime_error("ChemGillespieImpl::makeStep(): No Reaction was selected during the Gillespie step!");
        }
        
//        cout << "ChemGillespieImpl::makeStep(): mu=" << mu << ", rates_sum=" << rates_sum << ", _a_total=" << _a_total
//        << ", rn_selected*=" << rn_selected << endl;
//        rn_selected->printSelf();
        
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
        for(auto rit = r->dependents().begin(); rit!=r->dependents().end(); ++rit){
            RNodeGillespie *rn_other = static_cast<RNodeGillespie*>((*rit)->getRnode());
            rn_other->reComputePropensity();
            a_new = rn_other->getPropensity();
            a_penult = rn_other->getPenultStepPropensity();
            _a_total = _a_total - a_penult + a_new;
//            cout << "ChemGillespieImpl::_makeStep(): Dependents that were updated\n" << (*rn_other->getReaction());
        }
        
        // Send signals
#ifdef REACTION_SIGNALING
        r->emitSignal();
#endif
#ifdef RSPECIES_SIGNALING
//        for(auto sit = r->beginReactants(); sit!=r->endReactants(); ++sit){
//            if((*sit)->isSignaling())
//                (*sit)->emitSignal(-1);
//        }
//        for(auto sit = r->beginProducts(); sit!=r->endProducts(); ++sit){
//            if((*sit)->isSignaling())
//                (*sit)->emitSignal(1);
//        }
#endif
//        cout << "ChemGillespieImpl::_makeStep(): Ending..., _a_total=" << _a_total << "\n\n" << endl;
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
//        cout << "ChemGillespieImpl::activateReaction(ReactionBase *r): " << (*r);
        auto mit = _map_rnodes.find(r);
        if(mit!=_map_rnodes.end()){
            RNodeGillespie *rn_this = mit->second.get();
            rn_this->reComputePropensity();
            double a_new = rn_this->getPropensity();
            double a_penult = rn_this->getPenultStepPropensity();
            _a_total = _a_total - a_penult + a_new;
//            cout << "\nChemGillespieImpl::activateReaction(...): The following reaction has been activated: " << endl;
//            cout << (*r);
//            cout << "a_penult=" << a_penult << ", a_new=" << a_new << ", _a_total=" << _a_total << endl;
        }
        else
            throw std::out_of_range("ChemGillespieImpl::activateReaction(...): Reaction not found!");
    }
    
    void ChemGillespieImpl::passivateReaction(ReactionBase *r) {
//        cout << "ChemGillespieImpl::passivateReaction(ReactionBase *r): " << (*r);
        auto mit = _map_rnodes.find(r);
        if(mit==_map_rnodes.end())
            throw std::out_of_range("ChemGillespieImpl::passivateReaction(...): Reaction not found!");
        RNodeGillespie *rn_this = mit->second.get();
        
        double a_new, a_penult;
        a_penult = rn_this->getPropensity();
        rn_this->setPenultA(a_penult);
        a_new = 0;
        rn_this->setA(a_new);
        _a_total = _a_total - a_penult + a_new;

//        cout << "\nChemGillespieImpl::passivateReaction(...): The following reaction has been passivated: " << endl;
//        cout << (*r);
//        cout << "a_penult=" << a_penult << ", a_new=" << a_new << ", _a_total=" << _a_total << endl;
    }
    
} // end of namespace