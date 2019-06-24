
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

#include "ChemSimpleGillespieImpl.h"
#include "Rand.h"

#include "Rand.h"

void ChemSimpleGillespieImpl::initialize() {
    resetTime();
}


ChemSimpleGillespieImpl::~ChemSimpleGillespieImpl() noexcept {
    _reactions.clear();
}

floatingpoint ChemSimpleGillespieImpl::generateTau(floatingpoint a){

#if defined(_MSC_VER) && defined(_DEBUG) // MSVC requires the parameter to be positive
    if(a == 0.0) {
        _exp_distr.param(exponential_distribution<floatingpoint>::param_type(numeric_limits<floatingpoint>::min()));
        return numeric_limits<floatingpoint>::infinity();
    }
#endif

    exponential_distribution<floatingpoint>::param_type pm(a);
    _exp_distr.param(pm);
    return _exp_distr(Rand::eng);
}

floatingpoint ChemSimpleGillespieImpl::generateUniform(){
    return _uniform_distr(Rand::eng);
}

floatingpoint ChemSimpleGillespieImpl::computeTotalA(){
    floatingpoint rates_sum = 0;
    for (auto &r : _reactions){
        rates_sum+=r->computePropensity();
    }
    return rates_sum;
}

bool ChemSimpleGillespieImpl::makeStep() {
    
    floatingpoint a_total = computeTotalA();
    
    // this means that the network has come to a halt
    if(a_total<1e-15)
        return false;

    floatingpoint tau = generateTau(a_total);
    _t+=tau;
    syncGlobalTime();
    
    //Gillespie algorithm's second step: finding which reaction happened;
    floatingpoint mu = a_total*generateUniform();
    floatingpoint rates_sum = 0;
    ReactionBase* r_selected = nullptr;
    for (auto &r : _reactions){
        
        rates_sum+=r->computePropensity();
    
        if(rates_sum>mu){
            r_selected = r;
            break;
        }
    }
    if(r_selected==nullptr){
        cout << "ChemSimpleGillespieImpl::makeStep() for loop: rates_sum=" <<
            rates_sum << ", mu="
            << mu << ", a_total=" << a_total << endl;
        throw runtime_error( "ChemSimpleGillespieImpl::makeStep(): No Reaction was selected during the Gillespie step!");
    }
    r_selected->makeStep();

    // Send signal
#ifdef REACTION_SIGNALING
    r_selected->emitSignal();
#endif
    
    return true;
}

void ChemSimpleGillespieImpl::addReaction(ReactionBase *r) {
    auto vit = find(_reactions.begin(), _reactions.end(), r);
    if(vit==_reactions.end())
        _reactions.push_back(r);
}

void ChemSimpleGillespieImpl::removeReaction(ReactionBase *r) {
    auto vit = find(_reactions.begin(), _reactions.end(), r);
    if(vit!=_reactions.end())
        _reactions.erase(vit);
}

void ChemSimpleGillespieImpl::printReactions() const {
    for (auto &r : _reactions)
        cout << (*r);
}
