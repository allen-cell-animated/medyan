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
#include <boost/pool/pool.hpp>

//struct RSpeciesPoolTag { };

/// Print self into an iostream
std::ostream& operator<<(std::ostream& os, const chem::RSpecies& s){
    os << s.getFullName() << "[" << s.getN() << "]";
    return os;
}

namespace chem {
    
    using namespace std;
    boost::pool<> allocator_64bytes(sizeof(RSpecies),1024);
 
void* RSpecies::operator new(size_t size)
{
//    cout << "RSpecies::operator new(std::size_t size) called..." << endl;
    void *ptr = allocator_64bytes.malloc();
    return ptr;
    // RSpecies* cell = new (ptr) RSpecies();
}
    
void RSpecies::operator delete(void* ptr) noexcept
{
//    cout << "RSpecies::operator operator delete(void* ptr) called..." << endl;
    allocator_64bytes.free(ptr);
}

    
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
