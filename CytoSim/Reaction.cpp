//
//  Reaction.cpp
//  CytoSim-Experimenting
//
//  Created by Garegin Papoian on 4/19/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include <iostream>
#include "Reaction.h"

using namespace std;

std::vector<Reaction*> Reaction::_getAffectedReactions() {
    std::unordered_set<Reaction*> rxns;
    for(auto s : _species){
        rxns.insert(s->beginReactantReactions(),s->endReactantReactions());
    }
    //        std::sort(rxns.begin(),rxns.end());
    //        rxns.erase(std::unique(rxns.begin(),rxns.end()), rxns.end());
    rxns.erase(this);
    return std::vector<Reaction*>(rxns.begin(),rxns.end());
}

void Reaction::_registerNewDependent(Reaction *r){
    if(std::find(_dependents.begin(),_dependents.end(),r)==_dependents.end())
        _dependents.push_back(r);
}

void Reaction::_unregisterDependent(Reaction *r){
    auto it=std::find(_dependents.begin(),_dependents.end(),r);
    if(it!=_dependents.end())
        _dependents.erase(it);
}

Reaction::Reaction (std::initializer_list<Species*> species, unsigned char M, unsigned char N, float rate) : _species(species), _m(M), _rate(rate) {
    _species.shrink_to_fit();
    assert(_species.size()==(M+N) && "Reaction Ctor Bug");
    _dependents=_getAffectedReactions();
    for(auto s=_species.begin(); s<(_species.begin()+_m);++s)
    {
        for(auto r = (*s)->beginReactantReactions(); r!=(*s)->endReactantReactions(); ++r){
            if(this!=(*r))
                (*r)->_registerNewDependent(this);
        }
        for(auto r = (*s)->beginProductReactions(); r!=(*s)->endProductReactions(); ++r){
            if(this!=(*r))
                (*r)->_registerNewDependent(this);
        }
    }
    std::for_each(_species.begin(), _species.begin()+_m, [this](Species* s){s->addAsReactant(this);} );
    std::for_each(_species.begin()+_m, _species.end(),   [this](Species* s){s->addAsProduct(this);} );
}   

Reaction::~Reaction() {
    std::for_each(_species.begin(), _species.begin()+_m, [this](Species* s){s->removeAsReactant(this);} );
    std::for_each(_species.begin()+_m, _species.end(),   [this](Species* s){s->removeAsProduct(this);} );
    for(auto s=_species.begin(); s<(_species.begin()+_m);++s)
    {
        for(auto r = (*s)->beginReactantReactions(); r!=(*s)->endReactantReactions(); ++r){
            (*r)->_unregisterDependent(this);
//            cout << "Reaction::~Reaction(), unregistered dependent reaction, " << this << ", from, " << (*r) << endl; 
        }
        for(auto r = (*s)->beginProductReactions(); r!=(*s)->endProductReactions(); ++r){
            (*r)->_unregisterDependent(this);
//            cout << "Reaction::~Reaction(), unregistered dependent reaction, " << this << ", from, " << (*r) << endl; 
        }
    }
};

void Reaction::printSelf() {
    unsigned char i=0;
    for(auto s: _species){
        if(i==_m)
            cout << " ---> ";
        cout << s->getFullName() << "[" << (int)s->getN() << "]";
        if(i<_m-1)
            cout << " + ";
        ++i;
    }
    std::cout << ", " << "curr_rate = " << _rate << ", a=" << computePropensity() << ", Reaction ptr=" << this << "\n";
}
