//
//  ReactionBase.cpp
//  CytoSim
//
//  Created by Garegin Papoian on 9/18/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include "ReactionBase.h"
#include "Composite.h"

using namespace std;

namespace chem {
    
    ReactionBase::ReactionBase (float rate) :
    _rnode(nullptr), _rate(rate), _parent(nullptr), _rate_bare(rate)
    {
#ifdef REACTION_SIGNALING
        _signal=nullptr;
#endif
#if defined TRACK_ZERO_COPY_N || defined TRACK_UPPER_COPY_N
        _passivated=false;
#endif
    }
    
    Composite* ReactionBase::getRoot()
    {
        if(hasParent())
            return this->getParent()->getRoot();
        return nullptr;
    }
    
    
    void ReactionBase::registerNewDependent(ReactionBase *r){
        if(std::find(_dependents.begin(),_dependents.end(),r)==_dependents.end())
            _dependents.push_back(r);
    }
    
    void ReactionBase::unregisterDependent(ReactionBase *r){
        auto it=std::find(_dependents.begin(),_dependents.end(),r);
        //    cout << "ReactionBase::unregisterDependent: " << this << ", this rxn ptr needs to be erased from the dependent's list" << r << endl;
        if(it!=_dependents.end())
            _dependents.erase(it);
    }
    
#ifdef REACTION_SIGNALING
    void ReactionBase::startSignaling () {
        _signal = new ReactionEventSignal;
    }
    
    void ReactionBase::stopSignaling () {
        if (_signal!=nullptr)
            delete _signal;
        _signal = nullptr;
    }
    
    boost::signals2::connection ReactionBase::connect(std::function<void (ReactionBase *)> const &react_callback, int priority) {
        if (!isSignaling())
            startSignaling();
        return _signal->connect(priority, react_callback);
    }
#endif
    
    void ReactionBase::printDependents()  {
        cout << "ReactionBase: ptr=" << this << "\n"
        << (*this) << "the following ReactionBase objects are dependents: ";
        if(_dependents.size()==0)
            cout << "NONE" << endl;
        else
            cout << endl;
        for(auto r : _dependents)
            cout << (*r) << endl;
    }
}