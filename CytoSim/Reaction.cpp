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

Reaction::Reaction (std::initializer_list<Species*> species, unsigned char M, unsigned char N, float rate) : _species(species), _m(M), _rate(rate) {
    _species.shrink_to_fit();
    assert(_species.size()==(M+N) && "Reaction Ctor Bug");
    std::for_each(beginReactants(), endReactants(), [this](Species* s){s->addAsReactant(this);} );
    std::for_each(beginProducts(), endProducts(),   [this](Species* s){s->addAsProduct(this);} );
    _rnode=nullptr;
}   

Reaction::~Reaction() {
    std::for_each(beginReactants(), endReactants(), [this](Species* s){s->removeAsReactant(this);} );
    std::for_each(beginProducts(), endProducts(),   [this](Species* s){s->removeAsProduct(this);} );
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
