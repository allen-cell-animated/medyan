//
//  Reaction.cpp
//  CytoSim-Experimenting
//
//  Created by Garegin Papoian on 4/19/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include <boost/pool/pool.hpp>

#include <iostream>
#include "Reaction.h"
#include "ChemRNode.h"
#include "Composite.h"

using namespace std;

namespace chem {

boost::pool<> allocator_88bytes(sizeof(Reaction),1024);

void* Reaction::operator new(size_t size)
{
    //    cout << "Reaction::operator new(std::size_t size) called..." << endl;
    void *ptr = allocator_88bytes.malloc();
    return ptr;
    // RSpecies* cell = new (ptr) RSpecies();
}

void Reaction::operator delete(void* ptr) noexcept
{
    //    cout << "Reaction::operator operator delete(void* ptr) called..." << endl;
    allocator_88bytes.free(ptr);
}
    

    
Reaction::Reaction (std::initializer_list<Species*> species, unsigned char M, unsigned char N, float rate) : 
 _rnode(nullptr), _rate(rate), _parent(nullptr), _m(M)
    {
#ifdef REACTION_SIGNALING
    _signal=nullptr;
#endif
#if defined TRACK_ZERO_COPY_N || TRACK_UPPER_COPY_N
    _passivated=false;
#endif
    for( auto &s : species){
        _rspecies.push_back(&s->getRSpecies());
    }
    _rspecies.shrink_to_fit();
    assert(_rspecies.size()==(M+N) && "Reaction Ctor Bug");
    _dependents=getAffectedReactions();
//    cout << "Reaction::Reaction(...): " << this << endl;
//    for (auto rr : _dependents)
//        cout <<(*rr);
//    cout << endl;
//    activateReactionUnconditional();
    std::for_each(beginReactants(), endReactants(),
                  [this](RSpecies* s){s->addAsReactant(this);} );
    std::for_each(beginProducts(), endProducts(),   
                  [this](RSpecies* s){s->addAsProduct(this);} );
}
    
Reaction::Reaction (std::vector<Species*> species, unsigned char M, unsigned char N, float rate) :
_rnode(nullptr), _rate(rate), _parent(nullptr), _m(M)
    {
#ifdef REACTION_SIGNALING
    _signal=nullptr;
#endif
#if defined TRACK_ZERO_COPY_N || TRACK_UPPER_COPY_N
    _passivated=false;
#endif
    for( auto &s : species){
        _rspecies.push_back(&s->getRSpecies());
    }
    _rspecies.shrink_to_fit();
    assert(_rspecies.size()==(M+N) && "Reaction Ctor Bug");
    _dependents=getAffectedReactions();
    //    cout << "Reaction::Reaction(...): " << this << endl;
    //    for (auto rr : _dependents)
    //        cout <<(*rr);
    //    cout << endl;
    //    activateReactionUnconditional();
    std::for_each(beginReactants(), endReactants(),
                  [this](RSpecies* s){s->addAsReactant(this);} );
    std::for_each(beginProducts(), endProducts(),
                  [this](RSpecies* s){s->addAsProduct(this);} );
}

Reaction::~Reaction() noexcept
    {
    std::for_each(beginReactants(), endReactants(), [this](RSpecies* s){s->removeAsReactant(this);} );
    std::for_each(beginProducts(), endProducts(),   [this](RSpecies* s){s->removeAsProduct(this);} );
    // passivateReaction();
#ifdef REACTION_SIGNALING
    if(_signal!=nullptr)
        delete _signal;
#endif
    }

Composite* Reaction::getRoot()
{
    if(hasParent())
        return this->getParent()->getRoot();
    return nullptr;
}
    
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
    
void Reaction::unregisterDependent(Reaction *r){
    auto it=std::find(_dependents.begin(),_dependents.end(),r);
    //    cout << "Reaction::unregisterDependent: " << this << ", this rxn ptr needs to be erased from the dependent's list" << r << endl;
    if(it!=_dependents.end())
        _dependents.erase(it);
}

void Reaction::activateReactionUnconditional(){
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
#if defined TRACK_ZERO_COPY_N || TRACK_UPPER_COPY_N
    _passivated=false;
#endif
    if(_rnode!=nullptr)
        _rnode->activateReaction();
}

void Reaction::passivateReaction() {
    if(isPassivated())
        return;
    for(auto s=beginReactants(); s<endReactants();++s)
    {
        for(auto r = (*s)->beginReactantReactions(); r!=(*s)->endReactantReactions(); ++r){
            (*r)->unregisterDependent(this);
        }
        for(auto r = (*s)->beginProductReactions(); r!=(*s)->endProductReactions(); ++r){
            (*r)->unregisterDependent(this);
        }
    }
#if defined TRACK_ZERO_COPY_N || TRACK_UPPER_COPY_N
    _passivated=true;
#endif
    if(_rnode!=nullptr)
        _rnode->passivateReaction();
}

#ifdef REACTION_SIGNALING
void Reaction::startSignaling () {
    _signal = new ReactionEventSignal;
}

void Reaction::stopSignaling () {
    if (_signal!=nullptr)
        delete _signal;
    _signal = nullptr;
}

boost::signals2::connection Reaction::connect(std::function<void (Reaction *)> const &react_callback, int priority) {
    if (!isSignaling())
        startSignaling(); 
    return _signal->connect(priority, react_callback);
}
#endif

    
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
    cout << "Reaction: ptr=" << this << "\n"
         << (*this) << "the following Reaction objects are dependents: ";
    if(_dependents.size()==0)
        cout << "NONE" << endl;
    else
        cout << endl;
    for(auto r : _dependents)
        cout << (*r) << endl;
}
    
} // end of namespace
