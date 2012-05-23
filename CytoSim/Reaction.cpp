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

namespace chem {
    
Reaction::Reaction (std::initializer_list<Species*> species, unsigned char M, unsigned char N, float rate) : 
_rnode(nullptr), _rate(rate), _m(M), _is_signaling (false) {
    for( auto &s : species){
        _rspecies.push_back(&s->getRSpecies());
    }
    _rspecies.shrink_to_fit();
    assert(_rspecies.size()==(M+N) && "Reaction Ctor Bug");
    _dependents=getAffectedReactions();
    //activateReaction();
    std::for_each(beginReactants(), endReactants(), 
                  [this](RSpecies* s){s->addAsReactant(this);} );
    std::for_each(beginProducts(), endProducts(),   [this](RSpecies* s){s->addAsProduct(this);} );
}   

Reaction::~Reaction() {
    std::for_each(beginReactants(), endReactants(), [this](RSpecies* s){s->removeAsReactant(this);} );
    std::for_each(beginProducts(), endProducts(),   [this](RSpecies* s){s->removeAsProduct(this);} );
    passivateReaction();
};

std::vector<Reaction*> Reaction::getAffectedReactions() {
    std::unordered_set<Reaction*> rxns;
    for(auto s : _rspecies){
        rxns.insert(s->beginReactantReactions(),s->endReactantReactions());
    }
    //        std::sort(rxns.begin(),rxns.end());
    //        rxns.erase(std::unique(rxns.begin(),rxns.end()), rxns.end());
    rxns.erase(this);
    return std::vector<Reaction*>(rxns.begin(),rxns.end());
}
    
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


void Reaction::makeSignaling (ChemSignal &sm) {
    sm.addSignalingReaction(this);
    _is_signaling=true;
}

void Reaction::stopSignaling (ChemSignal &sm) {
    sm.disconnect(this);
    _is_signaling=false;
}

std::ostream& operator<<(std::ostream& os, const Reaction& rr){
    unsigned char i=0;
    for (auto sit = rr.cbeginReactants(); sit!=rr.cendReactants(); ++sit){
        os << (*sit)->getFullName() << "{" << (*sit)->getN() << "}";
        if(i<rr.getM()-1)
            os << " + ";
        ++i;
    }
    os << " ---> ";
    i=0;
    for (auto sit = rr.cbeginProducts(); sit!=rr.cendProducts(); ++sit){
        os << (*sit)->getFullName() << "{" << (*sit)->getN() << "}";
        if(i<((rr.getN()-rr.getM())-1))
            os << " + ";
        ++i;
    }
    os << ", " << "curr_rate = " << rr.getRate() << ", a=" <<rr.computePropensity() << ", Reaction ptr=" << &rr << "\n";
    return os;
}    

void Reaction::printDependents()  {
    cout << "Reaction: ptr=" << this << ", the following Reaction objects are dependents:\n";
    for(auto r : _dependents)
        cout << r << endl;
}
}
