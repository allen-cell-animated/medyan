//
//  RSpecies.cpp
//  CytoSim
//
//  Created by Garegin Papoian on 5/22/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include <iostream>
#include "RSpecies.h"
#include "Reaction.h"

#ifdef BOOST_MEM_POOL
    #include <boost/pool/pool.hpp>
    #include <boost/pool/pool_alloc.hpp>
#endif

//struct RSpeciesPoolTag { };

/// Print self into an iostream
std::ostream& operator<<(std::ostream& os, const chem::RSpecies& s){
    os << s.getFullName() << "[" << s.getN() << "]";
    return os;
}

namespace chem {
    
using namespace std;

#ifdef BOOST_MEM_POOL
boost::pool<> allocator_rspecies(sizeof(RSpecies),BOOL_POOL_NSIZE);
 
void* RSpecies::operator new(size_t size)
{
//    cout << "RSpecies::operator new(std::size_t size) called..." << endl;
//    void *ptr = allocator_rspecies.malloc();
    void *ptr = boost::fast_pool_allocator<RSpecies>::allocate();

    return ptr;
    // RSpecies* cell = new (ptr) RSpecies();
}
    
void RSpecies::operator delete(void* ptr) noexcept
{
//    cout << "RSpecies::operator operator delete(void* ptr) called..." << endl;
//    allocator_rspecies.free(ptr);
    boost::fast_pool_allocator<RSpecies>::deallocate((RSpecies*)ptr);
}
#endif
    
std::string RSpecies::getFullName() const {
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
    
}
