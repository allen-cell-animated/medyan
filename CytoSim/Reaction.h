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


#include "common.h"
#include "Species.h"

class RNode;

class Reaction {
private:
    std::vector<Species*> _species;
    std::vector<Reaction*> _dependents;
    RNode* _rnode;
    const unsigned char _m;
    float _rate;
    void _registerNewDependent(Reaction *r);
    void _unregisterDependent(Reaction *r);
public:
    // Constructor    
    Reaction (std::initializer_list<Species*> species, unsigned char M, unsigned char N, float rate);    
    Reaction (const Reaction &r) = delete; // no copying (including all derived classes)
    Reaction& operator=(Reaction &r) = delete;  // no assignment (including all derived classes)
    // Destructor
    ~Reaction();
//    // Setters 
    void setRate(float rate) {_rate=rate;}
    void setRnode(RNode *rhs) {_rnode=rhs;} 
//    // Accessors 
    float getRate() const {return _rate;}
    RNode* getRnode() const {return _rnode;} 
    //Iterators
    vr_iterator beginAffected() {return _dependents.begin();}
    vr_iterator endAffected() {return _dependents.end();}
    vsp_iterator beginReactants() {return _species.begin();}
    vsp_iterator endReactants() {return _species.begin()+_m;}
    vsp_iterator beginProducts() {return _species.begin()+_m;}
    vsp_iterator endProducts() {return _species.end();}

//    //Fire Reactions
    void makeStep() {
        for(auto sit = beginReactants(); sit!=endReactants(); ++sit) (*sit)->incrementN(-1);
        for(auto sit = beginProducts(); sit!=endProducts(); ++sit) (*sit)->incrementN(+1);
    }
    float computePropensity () {
        return std::accumulate(beginReactants(), endReactants(), 
                               _rate, 
                               [](float prod, Species *s){ 
                                   return prod*=s->getN();
                               } );
    }
    void printSelf () ; 
    void printDependents() ;
    std::vector<Reaction*> getAffectedReactions();
};


#endif
