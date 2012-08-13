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
    
    RNodeGillespie::RNodeGillespie(Reaction *r, ChemGillespieImpl &chem_Gillespie) :
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
        cout << "RNodeGillespie: ptr=" << this << ", a=" << _a << ", points to Reaction:\n";
        cout << (*_react);
    }
    
    void RNodeGillespie::printDependents() const {
        cout << "RNodeGillespie: ptr=" << this << ", the following RNodeGillespie objects are dependents:\n\n";
        for(auto rit = _react->dependents().begin(); rit!=_react->dependents().end(); ++rit){
            //        cout << "I am here [" << (*rit)->getRnode() << "]" << endl;
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
        _a_total = 0;

        for (auto &mit : _map_rnodes_inactive){
            _map_rnodes[mit.first]=std::move(mit.second);
        }
        _map_rnodes_inactive.clear();

        for (auto &x : _map_rnodes){
            auto rn = x.second.get();
            rn->getReaction()->activateReactionUnconditional();
            rn->reset();
            double a_new = rn->getPropensity();
            _a_total+=a_new;
        }
    }
    
    
    ChemGillespieImpl::~ChemGillespieImpl() {
        _map_rnodes.clear();
    }
    
    double ChemGillespieImpl::generateTau(double a){
        exponential_distribution<double>::param_type pm(a);
        _exp_distr.param(pm);
        return _exp_distr(_eng);
    }
    
    double ChemGillespieImpl::generateUniform(){
        return _uniform_distr(_eng);
    }
    
    bool ChemGillespieImpl::makeStep()
    {
//        cout << "\n\n[ChemGillespieImpl::_makeStep(): Starting..." << endl;
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
            rates_sum+=rn->getPropensity();
//            rn->printSelf();
//            cout << rn->getReaction()->computePropensity() << " " << rn->getPropensity() << endl;
//            if(abs(rn->getReaction()->computePropensity()-rn->getPropensity())>0.00001){
//                rn->printSelf();
//                cout << rn->getReaction()->computePropensity() << " " << rn->getPropensity() << endl;
//            }
//            cout << "ChemGillespieImpl::makeStep() for loop: a=" << rn->getPropensity() << ", rates_sum=" << rates_sum << ", mu="
//                << mu << ", _a_total=" << _a_total << endl;
            if(rates_sum>mu){
                rn_selected = rn;
                break;
            }
        }
        
//        cout << "ChemGillespieImpl::makeStep(): mu=" << mu << ", rates_sum=" << rates_sum << ", _a_total=" << _a_total
//        << ", rn_selected*=" << rn_selected << endl;
        
        if(rn_selected==nullptr){
            cout << "ChemGillespieImpl::makeStep() for loop: rates_sum=" << rates_sum << ", mu="
                        << mu << ", _a_total=" << _a_total << endl;
            throw std::runtime_error("ChemGillespieImpl::makeStep(): No Reaction was selected during the Gillespie step!");
        }
        
        rn_selected->makeStep();
        rn_selected->reComputePropensity();
        _a_total = _a_total - rn_selected->getPenultStepPropensity() + rn_selected->getPropensity();

//        if(!rn->isPassivated()){
//            rn->generateNewRandTau();
//            rn->updateHeap();
//        }
        
        //    cout << "ChemGillespieImpl::makeStep(): RNodeGillespie ptr=" << rn << " made a chemical step. t=" << _t << "\n" << endl;
        //    rn->printSelf();
        // Updating dependencies
        Reaction *r = rn_selected->getReaction();
        for(auto rit = r->dependents().begin(); rit!=r->dependents().end(); ++rit){
            RNodeGillespie *rn_other = static_cast<RNodeGillespie*>((*rit)->getRnode());
            rn_other->reComputePropensity();
//            double a_new = rn_other->getPropensity();
//            double a_penult = rn_other->getPenultStepPropensity();
            _a_total = _a_total - rn_other->getPenultStepPropensity() + rn_other->getPropensity();
        }
        
        // Send signals
        r->emitSignal();
        for(auto sit = r->beginReactants(); sit!=r->endReactants(); ++sit){
            if((*sit)->isSignaling())
                (*sit)->emitSignal(-1);
        }
        for(auto sit = r->beginProducts(); sit!=r->endProducts(); ++sit){
            if((*sit)->isSignaling())
                (*sit)->emitSignal(1);
        }
        //    cout << "ChemGillespieImpl::_makeStep(): Ending...]\n\n" << endl;
        syncGlobalTime();
        return true;
    }
    
    void ChemGillespieImpl::addReaction(Reaction *r) {
        _map_rnodes.emplace(r,make_unique<RNodeGillespie>(r,*this));
        ++_n_reacts;
    }
    
    void ChemGillespieImpl::removeReaction(Reaction *r) {
        _map_rnodes.erase(r);
        --_n_reacts;
    }
    
    void ChemGillespieImpl::printReactions() const {
        for (auto &x : _map_rnodes){
            auto rn = x.second.get();
            rn->printSelf();
        }
    }
    
    void ChemGillespieImpl::activateReaction(Reaction *r) {
        auto mit = _map_rnodes_inactive.find(r);
        if(mit==_map_rnodes_inactive.end())
            return;
        _map_rnodes[r]=std::move(mit->second);
        _map_rnodes_inactive.erase(r);
//        cout << "ChemGillespieImpl::activateReaction(...): The following reaction has been activated: " << endl;
//        cout << (*r) << endl;
    }
    
    void ChemGillespieImpl::passivateReaction(Reaction *r) {
        auto mit = _map_rnodes.find(r);
        if(mit==_map_rnodes.end())
            return;
        _map_rnodes_inactive[r]=std::move(mit->second);
        _map_rnodes.erase(r);
//        cout << "ChemGillespieImpl::passivateReaction(...): The following reaction has been passivated: " << endl;
//        cout << (*r) << endl;
    }
    
} // end of namespace