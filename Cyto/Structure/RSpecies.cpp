//
//  RSpecies.cpp
//  CytoSim
//
//  Created by Garegin Papoian on 5/22/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include "RSpecies.h"
#include "Reaction.h"

#ifdef BOOST_MEM_POOL
    #include <boost/pool/pool.hpp>
    #include <boost/pool/pool_alloc.hpp>
#endif

RSpecies::~RSpecies(){
    ///SHOULD BE UNCOMMENTED...
    if(!(_as_reactants.empty() and _as_products.empty())) {
        cout << "Major bug: RSpecies should not contain Reactions when being destroyed." << endl;
    }
    //assert((_as_reactants.empty() and _as_products.empty()) && "Major bug: RSpecies should not contain Reactions when being destroyed.");//Should not throw an exception from a destructor - that would be undefined behavior

#ifdef RSPECIES_SIGNALING
    if(_signal!=nullptr)
        delete _signal;
#endif
}


/// Print self into an iostream
ostream& operator<<(ostream& os, const RSpecies& s){
    os << s.getFullName() << "[" << s.getN() << "]";
    return os;
}

#ifdef BOOST_MEM_POOL
boost::pool<> allocator_rspecies(sizeof(RSpecies),BOOL_POOL_NSIZE);
 
void* RSpecies::operator new(size_t size) {
    void *ptr = boost::fast_pool_allocator<RSpecies>::allocate();
    return ptr;
}
    
void RSpecies::operator delete(void* ptr) noexcept {
    boost::fast_pool_allocator<RSpecies>::deallocate((RSpecies*)ptr);
}
#endif
    
string RSpecies::getFullName() const {
    return _species.getFullName();
}


void RSpecies::activateAssocReactantReactions() {
    for (auto &r : _as_reactants)
        r->activateReaction();
}
    
void RSpecies::activateAssocProductReactions() {
    for (auto &r : _as_products)
        r->activateReaction();
}

void RSpecies::passivateAssocReactantReactions() {
    for (auto &r : _as_reactants)
        r->passivateReaction();
}
    
void RSpecies::passivateAssocProductReactions() {
    for (auto &r : _as_products)
        r->passivateReaction();
}
    
    

#ifdef RSPECIES_SIGNALING    
void RSpecies::startSignaling () {
    _signal = new RSpeciesCopyNChangedSignal;
}

void RSpecies::stopSignaling () {
    if (_signal!=nullptr)
        delete _signal;
    _signal = nullptr;
}
#endif // of RSPECIES_SIGNALING
    
