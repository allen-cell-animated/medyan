//
//  Reaction.h
//  CytoSim-Experimenting
//
//  Created by Garegin Papoian on 4/19/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_Experimenting_Reaction_h
#define CytoSim_Experimenting_Reaction_h

#include <iostream>
#include <algorithm>
#include <utility>
#include <unordered_set>
#include <cmath>
#include <initializer_list>

#include <boost/flyweight.hpp>

#include "common.h"
#include "Species.h"


class ReactionBase {
protected:
    std::vector<Species*> _species;
    const unsigned char _m;
    float _rate;
public:
    ReactionBase (const ReactionBase &r) = delete; // no copying (including all derived classes)
    ReactionBase& operator=(ReactionBase &r) = delete;  // no assignment (including all derived classes)
//    // Setters 
    void setRate(float rate){_rate=rate;}
//    // Accessors 
    float getRate(){return _rate;}
//    //Fire Reactions
    void makeStep() {
        std::for_each(_species.begin(), _species.begin()+_m, [](Species* s){s->incrementN(-1);} );
        std::for_each(_species.begin()+_m, _species.end(),   [](Species* s){s->incrementN(+1);} );
    }
    float computePropensity (){
        return std::accumulate(_species.begin(), _species.begin()+_m, _rate, [](float prod, Species *s){ return prod*=s->getN();} );
    }
    void printSelf() {
        std::cout << "Reaction, ptr=" << this << ", current rate=" << _rate << "\n";
        this->printSelfImpl();
    }
    std::unordered_set<ReactionBase*> getAffectedReactions() {
        std::unordered_set<ReactionBase*> rxns;
        for(auto s : _species){
            rxns.insert(s->beginFReactions(),s->endFReactions());
        }
//        std::sort(rxns.begin(),rxns.end());
//        rxns.erase(std::unique(rxns.begin(),rxns.end()), rxns.end());
        rxns.erase(this);
        std::cout << "getAffectedReactions():" << rxns.size() << std::endl;
        return rxns;
    }
//    // (Implementation: Virtual Functions)
protected:
    virtual void printSelfImpl() {};
// Constructors
    ReactionBase (std::initializer_list<Species*> species, unsigned char M, unsigned char N, float rate) : _species(species), _m(M), _rate(rate) {
        _species.shrink_to_fit();
        assert(_species.size()==(M+N) && "ReactionBase Ctor Bug");
        std::for_each(_species.begin(), _species.begin()+_m, [this](Species* s){s->addFReaction(this);} );
        std::for_each(_species.begin()+_m, _species.end(),   [this](Species* s){s->addBReaction(this);} );
}    
//    //Destructor
    virtual ~ReactionBase() {
        std::for_each(_species.begin(), _species.begin()+_m, [this](Species* s){s->removeFReaction(this);} );
        std::for_each(_species.begin()+_m, _species.end(),   [this](Species* s){s->removeBReaction(this);} );
    };
};


template <unsigned char M,unsigned char N> 
class Reaction : public ReactionBase {
protected:
    void printSelfImpl() { // override
        for (auto &s : _species)
            s->printSelf();
    };
public:
    Reaction<M,N> (std::initializer_list<Species*> species, float rate): ReactionBase(species,M,N,rate) {}
};

#endif
