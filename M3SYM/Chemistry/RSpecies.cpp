
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "RSpecies.h"

#include "Reaction.h"

#ifdef BOOST_MEM_POOL
    #include <boost/pool/pool.hpp>
    #include <boost/pool/pool_alloc.hpp>
#endif

RSpecies::~RSpecies() noexcept{

    assert((_as_reactants.empty() and _as_products.empty())
    && "Major bug: RSpecies should not contain Reactions when being destroyed.");
#ifdef RSPECIES_SIGNALING
    if(_signal!=nullptr)
        delete _signal;
#endif
}

ostream& operator<<(ostream& os, const RSpecies& s){
    os << s.getFullName() << "[" << s.getN() << "]";
    return os;
}

#ifdef BOOST_MEM_POOL
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
    
