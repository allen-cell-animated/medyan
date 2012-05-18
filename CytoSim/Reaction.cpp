//
//  Reaction.cpp
//  CytoSim-Experimenting
//
//  Created by Garegin Papoian on 4/19/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include <iostream>
#include "Reaction.h"
#include "ChemNRMImpl.h"

using namespace std;

std::vector<Reaction*> Reaction::getAffectedReactions() {
    std::unordered_set<Reaction*> rxns;
    for(auto s : _species){
        rxns.insert(s->beginReactantReactions(),s->endReactantReactions());
    }
    //        std::sort(rxns.begin(),rxns.end());
    //        rxns.erase(std::unique(rxns.begin(),rxns.end()), rxns.end());
    rxns.erase(this);
    return std::vector<Reaction*>(rxns.begin(),rxns.end());
}


Reaction::Reaction (std::initializer_list<Species*> species, unsigned char M, unsigned char N, float rate, bool is_signaling) : 
_species(species), _rate(rate), _m(M), _is_signaling (is_signaling) {
    _species.shrink_to_fit();
    assert(_species.size()==(M+N) && "Reaction Ctor Bug");
    _dependents=getAffectedReactions();
    //activateReaction();
    std::for_each(beginReactants(), endReactants(), [this](Species* s){s->addAsReactant(this);} );
    std::for_each(beginProducts(), endProducts(),   [this](Species* s){s->addAsProduct(this);} );
    _rnode=nullptr;
}   

Reaction::~Reaction() {
    std::for_each(beginReactants(), endReactants(), [this](Species* s){s->removeAsReactant(this);} );
    std::for_each(beginProducts(), endProducts(),   [this](Species* s){s->removeAsProduct(this);} );
    passivateReaction();
};

void Reaction::registerNewDependent(Reaction *r){
    if(std::find(_dependents.begin(),_dependents.end(),r)==_dependents.end())
        _dependents.push_back(r);
}

void Reaction::activateReaction(){
    if(getProductOfReactants()==0) // One of the reactants is still at zero copy n, no need to activate yet...
        return;
    for(auto s=beginReactants(); s<endReactants();++s)
    {
        for(auto r = (*s)->beginReactantReactions(); r!=(*s)->endReactantReactions(); ++r){
            if(this!=(*r))
                (*r)->registerNewDependent(this);
        }
        for(auto r = (*s)->beginProductReactions(); r!=(*s)->endProductReactions(); ++r){
            if(this!=(*r))
                (*r)->registerNewDependent(this);
        }
    }
    if(_rnode!=nullptr)
        _rnode->activateReaction();
}

void Reaction::unregisterDependent(Reaction *r){
    auto it=std::find(_dependents.begin(),_dependents.end(),r);
    if(it!=_dependents.end())
        _dependents.erase(it);
}

void Reaction::passivateReaction() {
    for(auto s=beginReactants(); s<endReactants();++s)
    {
        for(auto r = (*s)->beginReactantReactions(); r!=(*s)->endReactantReactions(); ++r){
            (*r)->unregisterDependent(this);
        }
        for(auto r = (*s)->beginProductReactions(); r!=(*s)->endProductReactions(); ++r){
            (*r)->unregisterDependent(this);
        }
    }    
    if(_rnode!=nullptr)
        _rnode->passivateReaction();
}


void Reaction::makeSignaling (SignalingManager &sm) {
    sm.addSignalingReaction(this);
    _is_signaling=true;
}

void Reaction::stopSignaling (SignalingManager &sm) {
    sm.disconnect_semiprivate(this);
    _is_signaling=false;
}

void Reaction::printSelf() {
    unsigned char i=0;
    for(auto s: _species){
        if(i==_m)
            cout << " ---> ";
        cout << s->getFullName() << "{" << (int)s->getN() << "}";
        if(i<_m-1 || (i>=_m && i<_species.size()-1))
            cout << " + ";
        ++i;
    }
    std::cout << ", " << "curr_rate = " << _rate << ", a=" << computePropensity() << ", Reaction ptr=" << this << "\n";
}

void Reaction::printDependents()  {
    cout << "Reaction: ptr=" << this << ", the following Reaction objects are dependents:\n";
    for(auto r : _dependents)
        r->printSelf();
}

