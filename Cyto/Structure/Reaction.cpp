//
//  Reaction.cpp
//  CytoSim-Experimenting
//
//  Created by Garegin Papoian on 4/19/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include <iostream>
#include "Reaction.h"
#include "ChemRNode.h"
#include "SpeciesContainer.h"

#ifdef BOOST_MEM_POOL
#include <boost/pool/pool.hpp>
#include <boost/pool/pool_alloc.hpp>
#endif

using namespace std;

template <unsigned short M, unsigned short N>
    void Reaction<M,N>::activateReactionUnconditionalImpl(){
    for(auto i=0U; i<M; ++i)
    {
        RSpecies *s = _rspecies[i];
        for(auto r = s->beginReactantReactions(); r!=s->endReactantReactions(); ++r){
            if(this!=(*r))
                (*r)->registerNewDependent(this);
        }
        for(auto r = s->beginProductReactions(); r!=s->endProductReactions(); ++r){
            if(this!=(*r))
                (*r)->registerNewDependent(this);
        }
    }
#if defined TRACK_ZERO_COPY_N || defined TRACK_UPPER_COPY_N
    _passivated=false;
#endif
    if(_rnode!=nullptr)
        _rnode->activateReaction();
}

template <unsigned short M, unsigned short N>
void Reaction<M,N>::passivateReactionImpl() {
    if(isPassivated())
        return;
    for(auto i=0U; i<M; ++i)
    {
        RSpecies *s = _rspecies[i];
        for(auto r = s->beginReactantReactions(); r!=s->endReactantReactions(); ++r){
            (*r)->unregisterDependent(this);
        }
        for(auto r = s->beginProductReactions(); r!=s->endProductReactions(); ++r){
            (*r)->unregisterDependent(this);
        }
    }
#if defined TRACK_ZERO_COPY_N || defined TRACK_UPPER_COPY_N
    _passivated=true;
#endif
    if(_rnode!=nullptr)
        _rnode->passivateReaction();
}
  
template <unsigned short M, unsigned short N>
Reaction<M,N>* Reaction<M,N>::cloneImpl(const SpeciesPtrContainerVector &spcv)
{
    std::vector<Species*> species;
    for(auto &rs : _rspecies){
        int molec = rs->getSpecies().getMolecule();
        
        ///If species bulk, just add it
        Species* s = &rs->getSpecies();
        if(dynamic_cast<SpeciesBulk*>(s) != nullptr)
            species.push_back(s);
        
        ///otherwise, check if that species exists in the compartment
        else {
            auto vit = std::find_if(spcv.species().cbegin(),spcv.species().cend(),
                                    [molec](const std::unique_ptr<Species> &us){return us->getMolecule()==molec;});
            if(vit==spcv.species().cend())
                throw std::runtime_error("ReactionBase::Clone(): Species is not present.");
            species.push_back(vit->get());
        }
    }
    ///Create new reaction, copy ownership of signal
    Reaction* newReaction = new Reaction<M,N>(species,_rate);
    newReaction->_signal = _signal;

    return newReaction;
}
    
#ifdef BOOST_MEM_POOL
// boost::pool<> allocator_reaction(sizeof(Reaction<1,1>),BOOL_POOL_NSIZE);
template <unsigned short M, unsigned short N>
void* Reaction<M,N>::operator new(size_t size)
{
    //    cout << "Reaction<M,N>::operator new(std::size_t size) called..." << endl;
    //    void *ptr = allocator_ReactionBase.malloc();
    void *ptr = boost::fast_pool_allocator<Reaction<M,N>>::allocate();
    return ptr;
    // RSpecies* cell = new (ptr) RSpecies();
}

template <unsigned short M, unsigned short N>
void Reaction<M,N>::operator delete(void* ptr) noexcept
{
    //    cout << "ReactionBase::operator operator delete(void* ptr) called..." << endl;
    //    allocator_ReactionBase.free(ptr);
    boost::fast_pool_allocator<Reaction<M,N>>::deallocate((Reaction<M,N>*)ptr);
}
#endif
    
#ifdef BOOST_MEM_POOL
template void* Reaction<1,1>::operator new(size_t size);
template void Reaction<1,1>::operator delete(void* ptr);
template void* Reaction<2,1>::operator new(size_t size);
template void Reaction<2,1>::operator delete(void* ptr);
template void* Reaction<1,2>::operator new(size_t size);
template void Reaction<1,2>::operator delete(void* ptr);
template void* Reaction<2,2>::operator new(size_t size);
template void Reaction<2,2>::operator delete(void* ptr);
template void* Reaction<2,0>::operator new(size_t size);
template void Reaction<2,0>::operator delete(void* ptr);
template void* Reaction<1,3>::operator new(size_t size);
template void Reaction<1,3>::operator delete(void* ptr);
template void* Reaction<2,3>::operator new(size_t size);
template void Reaction<2,3>::operator delete(void* ptr);
template void* Reaction<3,2>::operator new(size_t size);
template void Reaction<3,2>::operator delete(void* ptr);
#endif

template void Reaction<1,1>::activateReactionUnconditionalImpl();
template void Reaction<1,1>::passivateReactionImpl();
template Reaction<1,1>* Reaction<1,1>::cloneImpl(const SpeciesPtrContainerVector &spcv);

template void Reaction<2,1>::activateReactionUnconditionalImpl();
template void Reaction<2,1>::passivateReactionImpl();
template Reaction<2,1>* Reaction<2,1>::cloneImpl(const SpeciesPtrContainerVector &spcv);

template void Reaction<1,2>::activateReactionUnconditionalImpl();
template void Reaction<1,2>::passivateReactionImpl();
template Reaction<1,2>* Reaction<1,2>::cloneImpl(const SpeciesPtrContainerVector &spcv);

template void Reaction<2,2>::activateReactionUnconditionalImpl();
template void Reaction<2,2>::passivateReactionImpl();
template Reaction<2,2>* Reaction<2,2>::cloneImpl(const SpeciesPtrContainerVector &spcv);

template void Reaction<2,0>::activateReactionUnconditionalImpl();
template void Reaction<2,0>::passivateReactionImpl();
template Reaction<2,0>* Reaction<2,0>::cloneImpl(const SpeciesPtrContainerVector &spcv);

template void Reaction<1,3>::activateReactionUnconditionalImpl();
template void Reaction<1,3>::passivateReactionImpl();
template Reaction<1,3>* Reaction<1,3>::cloneImpl(const SpeciesPtrContainerVector &spcv);

template void Reaction<2,3>::activateReactionUnconditionalImpl();
template void Reaction<2,3>::passivateReactionImpl();
template Reaction<2,3>* Reaction<2,3>::cloneImpl(const SpeciesPtrContainerVector &spcv);

template void Reaction<3,2>::activateReactionUnconditionalImpl();
template void Reaction<3,2>::passivateReactionImpl();
template Reaction<3,2>* Reaction<3,2>::cloneImpl(const SpeciesPtrContainerVector &spcv);

template void Reaction<3,0>::activateReactionUnconditionalImpl();
template void Reaction<3,0>::passivateReactionImpl();
template Reaction<3,0>* Reaction<3,0>::cloneImpl(const SpeciesPtrContainerVector &spcv);

