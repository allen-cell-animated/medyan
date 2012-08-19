//
//  ChemSimpleGillespieImpl.cpp
//  CytoSim
//
//  Created by Garegin Papoian on 8/16/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include "ChemSimpleGillespieImpl.h"

#include <iostream>

using namespace std;

namespace chem {
    
    void ChemSimpleGillespieImpl::initialize() {
        resetTime();
    }
    
    
    ChemSimpleGillespieImpl::~ChemSimpleGillespieImpl() {
        _reactions.clear();
    }
    
    double ChemSimpleGillespieImpl::generateTau(double a){
        exponential_distribution<double>::param_type pm(a);
        _exp_distr.param(pm);
        return _exp_distr(_eng);
    }
    
    double ChemSimpleGillespieImpl::generateUniform(){
        return _uniform_distr(_eng);
    }
    
    double ChemSimpleGillespieImpl::computeTotalA(){
        double rates_sum = 0;
        for (auto &r : _reactions){
            rates_sum+=r->computePropensity();
        }
        return rates_sum;
    }
    
    bool ChemSimpleGillespieImpl::makeStep()
    {
//        cout << "\n\n[ChemSimpleGillespieImpl::_makeStep(): Starting..." << endl;
//        cout << "Size of the reaction network, "<< _reactions.size() << endl;
//        printReactions();
        
        double a_total = computeTotalA();

        double tau = generateTau(a_total);
        _t+=tau;
        syncGlobalTime();
        
        //Gillespie algorithm's second step: finding which reaction happened;
        double mu = a_total*generateUniform();
        double rates_sum = 0;
        Reaction* r_selected = nullptr;
        for (auto &r : _reactions){
            rates_sum+=r->computePropensity();
            if(rates_sum>mu){
                r_selected = r;
                break;
            }
        }
        
        if(r_selected==nullptr){
            cout << "ChemSimpleGillespieImpl::makeStep() for loop: rates_sum=" << rates_sum << ", mu="
            << mu << ", a_total=" << a_total << endl;
            throw std::runtime_error("ChemSimpleGillespieImpl::makeStep(): No Reaction was selected during the Gillespie step!");
        }
        
//        cout << "ChemSimpleGillespieImpl::makeStep(): mu=" << mu << ", rates_sum=" << rates_sum << ", a_total=" << a_total
//        << ", r_selected*=" << (*r_selected) << endl;
        
        r_selected->makeStep();

        // Send signals
        r_selected->emitSignal();
        for(auto sit = r_selected->beginReactants(); sit!=r_selected->endReactants(); ++sit){
            if((*sit)->isSignaling())
                (*sit)->emitSignal(-1);
        }
        for(auto sit = r_selected->beginProducts(); sit!=r_selected->endProducts(); ++sit){
            if((*sit)->isSignaling())
                (*sit)->emitSignal(1);
        }
//        cout << "ChemSimpleGillespieImpl::_makeStep(): Ending..., _a_total=" << a_total << "\n\n" << endl;
        syncGlobalTime();
        return true;
    }
    
    void ChemSimpleGillespieImpl::addReaction(Reaction *r) {
        std::vector<Reaction*>::iterator vit = std::find(_reactions.begin(), _reactions.end(), r);
        if(vit==_reactions.end())
            _reactions.push_back(r);
    }
    
    void ChemSimpleGillespieImpl::removeReaction(Reaction *r) {
        std::vector<Reaction*>::iterator vit = std::find(_reactions.begin(), _reactions.end(), r);
        if(vit!=_reactions.end())
            _reactions.erase(vit);
    }
    
    void ChemSimpleGillespieImpl::printReactions() const {
        for (auto &r : _reactions){
            cout << (*r); 
        }
    }
        
} // end of namespace