
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "ReactionBase.h"

#include "Composite.h"

#include "ChemNRMImpl.h"

size_t ReactionBase::_Idcounter = 0;

ReactionBase::ReactionBase (float rate, bool isProtoCompartment, floatingpoint volumeFrac, int rateVolumeDepExp)
    : _rnode(nullptr), _parent(nullptr), _rate(rate),
    _rate_bare(rate), _isProtoCompartment(isProtoCompartment),
    _volumeFrac(volumeFrac), _rateVolumeDepExp(rateVolumeDepExp) {

	for(uint i = 0; i < RateMulFactorType::RATEMULFACTSIZE; i++)
		_ratemulfactors[i] = 1.0;

    // Scale the rate
	recalcRateVolumeFactor();
#ifdef REACTION_SIGNALING
    _signal=nullptr;
#endif
    //All reactions are generated passivated.
#if defined TRACK_ZERO_COPY_N || defined TRACK_UPPER_COPY_N
    _passivated=true;
#endif
    _Id = _Idcounter;
    _Idcounter++;
}

Composite* ReactionBase::getRoot() {
    if(hasParent())
        return this->getParent()->getRoot();
    return nullptr;
}

void ReactionBase::registerNewDependent(ReactionBase *r){ _dependents.insert(r);}

void ReactionBase::unregisterDependent(ReactionBase *r){ _dependents.erase(r);}

#ifdef REACTION_SIGNALING
void ReactionBase::startSignaling () {
    _signal = unique_ptr<ReactionEventSignal>(new ReactionEventSignal);
}

void ReactionBase::stopSignaling () {
    _signal = nullptr;
}

boost::signals2::connection ReactionBase::connect(
    function<void (ReactionBase *)> const &react_callback, int priority) {
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

bool afterchemsiminit = false;
void ReactionBase::activateReaction() {
#ifdef TRACK_ZERO_COPY_N
	if(areEqual(getProductOfReactants(), 0.0)) // One of the reactants is still at zero copy n,
		// no need to activate yet...
		return;
#endif
#ifdef TRACK_UPPER_COPY_N
	if(areEqual(getProductOfProducts(), 0.0)) // One of the products is at the maximum allowed
		//copy number, no need to activate yet...
		return;
#endif
	activateReactionUnconditional();
}

