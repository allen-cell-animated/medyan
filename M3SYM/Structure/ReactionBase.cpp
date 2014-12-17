
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "ReactionBase.h"

#include "Composite.h"

ReactionBase::ReactionBase (float rate, bool isProtoCompartment) :
_rnode(nullptr), _rate(rate), _parent(nullptr), _rate_bare(rate), _isProtoCompartment(isProtoCompartment)
{
#ifdef REACTION_SIGNALING
    _signal=nullptr;
#endif
#if defined TRACK_ZERO_COPY_N || defined TRACK_UPPER_COPY_N
    _passivated=true;
#endif
}

Composite* ReactionBase::getRoot() {
    if(hasParent())
        return this->getParent()->getRoot();
    return nullptr;
}


void ReactionBase::registerNewDependent(ReactionBase *r){
    if(find(_dependents.begin(),_dependents.end(),r)==_dependents.end())
        _dependents.push_back(r);
}

void ReactionBase::unregisterDependent(ReactionBase *r){
    auto it=find(_dependents.begin(),_dependents.end(),r);
    if(it!=_dependents.end())
        _dependents.erase(it);
}

#ifdef REACTION_SIGNALING
void ReactionBase::startSignaling () {
    _signal = unique_ptr<ReactionEventSignal>(new ReactionEventSignal);
}

void ReactionBase::stopSignaling () {
    _signal = nullptr;
}

boost::signals2::connection ReactionBase::connect(function<void (ReactionBase *)> const &react_callback, int priority) {
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

