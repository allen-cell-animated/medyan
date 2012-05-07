//
//  ChemNRMImpl.cpp
//  CytoSim
//
//  Created by Garegin Papoian on 5/6/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include <iostream>

using namespace std;

#include "ChemNRMImpl.h"



ReactionNodeNRM::ReactionNodeNRM(Reaction &r,  boost_heap &bh) : _dependents(0), _heap(bh), _react (r) {
    std::cout << "ReactionNodeNRM ctor, ptr=" << this << ", self handle, " << std::endl;
    _react.setRnode(this);
    auto rdeps=_react.getAffectedReactions();
    for(auto r : rdeps){
        auto rn = r->getRnode();
        if(rn!=nullptr)
            _dependents.push_back(rn);
    }
    std::cout << "Deps: " << _dependents.size() << std::endl; 
    for(auto s = _react.beginReactants(); s!=_react.endReactants(); ++s){
        for(auto r = (*s)->beginReactantReactions(); r!=(*s)->endReactantReactions(); ++r){
            auto rn = (*r)->getRnode();
            if(this!=rn and rn!=nullptr)
                rn->registerNewDependent(this);
        }
        for(auto r = (*s)->beginProductReactions(); r!=(*s)->endProductReactions(); ++r){
            auto rn = (*r)->getRnode();
            if(this!=rn and rn!=nullptr)
                rn->registerNewDependent(this);
        }
    }
    PQNode pqn = {this};
    pqn.tau=std::numeric_limits<float>::quiet_NaN();
    _handle = _heap.push(pqn);
}

ReactionNodeNRM::~ReactionNodeNRM(){
    for(auto s = _react.beginReactants(); s!=_react.endReactants(); ++s){
        {
            for(auto r = (*s)->beginReactantReactions(); r!=(*s)->endReactantReactions(); ++r){
                auto rn = (*r)->getRnode();
                if(rn != nullptr)
                    rn->unregisterDependent(this);
//                cout << "Reaction::~Reaction(), unregistered dependent reaction, " << this << ", from, " << rn << endl; 
            }
            for(auto r = (*s)->beginProductReactions(); r!=(*s)->endProductReactions(); ++r){
                auto rn = (*r)->getRnode();
                if(rn != nullptr)
                    rn->unregisterDependent(this);
//                cout << "Reaction::~Reaction(), unregistered dependent reaction, " << this << ", from, " << rn << endl; 
            }
        }
        _react.setRnode(nullptr);
    }
    _heap.erase(_handle);
}

void ReactionNodeNRM::makeStep(float t) {
    _react.makeStep();
    // auto a=_react.computePropensity();
    //randDrawTau();
    //updateNode();
    for(auto rit=_react.beginAffected();rit!=_react.endAffected();++rit){
        // get a_prev, tau_prev from the loop RNode
        // ...
//        float a_prev=1.0; //fake
//        float tau_prev=1.0;//fake
//        float a_new = (*rit)->computePropensity();
//        float tau_new = (a_prev/a_new)*(tau_prev-t)+t;
        // set tau_new ,a_new to the loop RNode
        // updateNode() the loop RNode
    }
}

void ReactionNodeNRM::updateHeap(){_heap.update(_handle);}
float ReactionNodeNRM::getTau() const {return (*_handle).tau;}
handle_t& ReactionNodeNRM::getHandle() {return _handle;}
void ReactionNodeNRM::setTau(float tau){(*_handle).tau=tau;}


void ReactionNodeNRM::registerNewDependent(ReactionNodeNRM *rn){
    auto it = _dependents.end();
    if(std::find(_dependents.begin(),_dependents.end(),rn)==_dependents.end()){
        _dependents.push_back(rn);
    }
}

void ReactionNodeNRM::unregisterDependent(ReactionNodeNRM *rn){
    auto it=std::find(_dependents.begin(),_dependents.end(),rn);
    if(it!=_dependents.end()){
        _dependents.erase(it);
    }
}

void ReactionNodeNRM::printSelf() const{
    std::cout << "ReactionNodeNRM: ptr=" << this << ", Reaction:\n";
    _react.printSelf();
}

void ReactionNodeNRM::printDependents() const {
    std::cout << "The following ReactionNodeNRM objects are dependents\n";
    for(auto r : _dependents)
        r->getReaction().printSelf();
}

//
//template <typename PriorityQueue>
//class ChemNRM {
//    ChemNRM(){_init();}
//private:
//    void _init();
//    const PriorityQueue* _pq;
//};
