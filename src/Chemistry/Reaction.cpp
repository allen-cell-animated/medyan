
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.0
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "Reaction.h"

#include "ChemRNode.h"
#include "SpeciesContainer.h"

#ifdef BOOST_MEM_POOL
#include <boost/pool/pool.hpp>
#include <boost/pool/pool_alloc.hpp>
#endif

template<unsigned short M, unsigned short N>
    void Reaction<M,N>::updatePropensityImpl() {

    //just update the rnode if not passivated
    if(_rnode!=nullptr && !_passivated) _rnode->activateReaction();
}


template <unsigned short M, unsigned short N>
    void Reaction<M,N>::activateReactionUnconditionalImpl(){
#ifdef TRACK_DEPENDENTS
    for(auto i=0U; i<M; ++i)
    {
        RSpecies *s = _rspecies[i];
        for(auto r = s->beginReactantReactions();
                 r!= s->endReactantReactions(); ++r){
            if(this!=(*r)) (*r)->registerNewDependent(this);
        }
        for(auto r = s->beginProductReactions();
                 r!= s->endProductReactions(); ++r){
            if(this!=(*r)) (*r)->registerNewDependent(this);
        }
    }
#endif
#if defined TRACK_ZERO_COPY_N || defined TRACK_UPPER_COPY_N
    _passivated=false;
#endif
    if(_rnode!=nullptr) _rnode->activateReaction();
}

template <unsigned short M, unsigned short N>
void Reaction<M,N>::passivateReactionImpl() {
    if(isPassivated()) return;
#ifdef TRACK_DEPENDENTS
    for(auto i=0U; i<M; ++i)
    {
        RSpecies *s = _rspecies[i];
        for(auto r = s->beginReactantReactions();
                 r!=s->endReactantReactions(); ++r){
            (*r)->unregisterDependent(this);
        }
        for(auto r = s->beginProductReactions();
                 r!=s->endProductReactions(); ++r){
            (*r)->unregisterDependent(this);
        }
    }
#endif
#if defined TRACK_ZERO_COPY_N || defined TRACK_UPPER_COPY_N
    _passivated=true;
#endif
    if(_rnode!=nullptr) _rnode->passivateReaction();
}
  
template <unsigned short M, unsigned short N>
Reaction<M,N>* Reaction<M,N>::cloneImpl(const SpeciesPtrContainerVector &spcv)
{
    vector<Species*> species;
    for(auto &rs : _rspecies){
        int molec = rs->getSpecies().getMolecule();
        
        //check if that species exists in the compartment
        auto vit = find_if(spcv.species().cbegin(),spcv.species().cend(),
           [molec](const unique_ptr<Species> &us){return us->getMolecule()==molec;});
        
        //if we didn't find it, use the old species
        if(vit==spcv.species().cend()) species.push_back(&rs->getSpecies());
        else species.push_back(vit->get());
    }
    //Create new reaction, copy ownership of signal
    Reaction* newReaction = new Reaction<M,N>(species, _rate_bare);
    newReaction->setRate(_rate);
    
#ifdef REACTION_SIGNALING
    newReaction->_signal = std::move(_signal);
    _signal = nullptr;
#endif
    //Copy reaction type
    newReaction->_reactionType = _reactionType;
    return newReaction;
}

void DiffusionReaction::updatePropensityImpl() {
    
    //just update the rnode if not passivated
    if(_rnode!=nullptr && !_passivated) _rnode->activateReaction();
}

    
#ifdef BOOST_MEM_POOL
template <unsigned short M, unsigned short N>
void* Reaction<M,N>::operator new(size_t size) {
    void *ptr = boost::fast_pool_allocator<Reaction<M,N>>::allocate();
    return ptr;
}

template <unsigned short M, unsigned short N>
void Reaction<M,N>::operator delete(void* ptr) noexcept {
    boost::fast_pool_allocator<Reaction<M,N>>::deallocate((Reaction<M,N>*)ptr);
}

void* DiffusionReaction::operator new(size_t size) {
    void *ptr = boost::fast_pool_allocator<DiffusionReaction>::allocate();
    return ptr;
}

void DiffusionReaction::operator delete(void* ptr) noexcept {
     boost::fast_pool_allocator<DiffusionReaction>::deallocate((DiffusionReaction*)ptr);
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
template void* Reaction<3,1>::operator new(size_t size);
template void Reaction<3,1>::operator delete(void* ptr);
template void* Reaction<3,0>::operator new(size_t size);
template void Reaction<3,0>::operator delete(void* ptr);
template void* Reaction<2,4>::operator new(size_t size);
template void Reaction<2,4>::operator delete(void* ptr);
template void* Reaction<2,5>::operator new(size_t size);
template void Reaction<2,5>::operator delete(void* ptr);
template void* Reaction<4,0>::operator new(size_t size);
template void Reaction<4,0>::operator delete(void* ptr);
template void* Reaction<5,2>::operator new(size_t size);
template void Reaction<5,2>::operator delete(void* ptr);
template void* Reaction<4,2>::operator new(size_t size);
template void Reaction<4,2>::operator delete(void* ptr);
#endif

template void Reaction<1,1>::updatePropensityImpl();
template void Reaction<1,1>::activateReactionUnconditionalImpl();
template void Reaction<1,1>::passivateReactionImpl();
template Reaction<1,1>* Reaction<1,1>::cloneImpl(const SpeciesPtrContainerVector &spcv);

template void Reaction<2,1>::updatePropensityImpl();
template void Reaction<2,1>::activateReactionUnconditionalImpl();
template void Reaction<2,1>::passivateReactionImpl();
template Reaction<2,1>* Reaction<2,1>::cloneImpl(const SpeciesPtrContainerVector &spcv);

template void Reaction<1,2>::updatePropensityImpl();
template void Reaction<1,2>::activateReactionUnconditionalImpl();
template void Reaction<1,2>::passivateReactionImpl();
template Reaction<1,2>* Reaction<1,2>::cloneImpl(const SpeciesPtrContainerVector &spcv);

template void Reaction<2,2>::updatePropensityImpl();
template void Reaction<2,2>::activateReactionUnconditionalImpl();
template void Reaction<2,2>::passivateReactionImpl();
template Reaction<2,2>* Reaction<2,2>::cloneImpl(const SpeciesPtrContainerVector &spcv);

template void Reaction<2,0>::updatePropensityImpl();
template void Reaction<2,0>::activateReactionUnconditionalImpl();
template void Reaction<2,0>::passivateReactionImpl();
template Reaction<2,0>* Reaction<2,0>::cloneImpl(const SpeciesPtrContainerVector &spcv);

template void Reaction<1,3>::updatePropensityImpl();
template void Reaction<1,3>::activateReactionUnconditionalImpl();
template void Reaction<1,3>::passivateReactionImpl();
template Reaction<1,3>* Reaction<1,3>::cloneImpl(const SpeciesPtrContainerVector &spcv);

template void Reaction<2,3>::updatePropensityImpl();
template void Reaction<2,3>::activateReactionUnconditionalImpl();
template void Reaction<2,3>::passivateReactionImpl();
template Reaction<2,3>* Reaction<2,3>::cloneImpl(const SpeciesPtrContainerVector &spcv);

template void Reaction<3,2>::activateReactionUnconditionalImpl();
template void Reaction<3,2>::passivateReactionImpl();
template Reaction<3,2>* Reaction<3,2>::cloneImpl(const SpeciesPtrContainerVector &spcv);

template void Reaction<3,1>::updatePropensityImpl();
template void Reaction<3,1>::activateReactionUnconditionalImpl();
template void Reaction<3,1>::passivateReactionImpl();
template Reaction<3,1>* Reaction<3,1>::cloneImpl(const SpeciesPtrContainerVector &spcv);

template void Reaction<3,0>::updatePropensityImpl();
template void Reaction<3,0>::activateReactionUnconditionalImpl();
template void Reaction<3,0>::passivateReactionImpl();
template Reaction<3,0>* Reaction<3,0>::cloneImpl(const SpeciesPtrContainerVector &spcv);

template void Reaction<2,4>::updatePropensityImpl();
template void Reaction<2,4>::activateReactionUnconditionalImpl();
template void Reaction<2,4>::passivateReactionImpl();
template Reaction<2,4>* Reaction<2,4>::cloneImpl(const SpeciesPtrContainerVector &spcv);

template void Reaction<2,5>::updatePropensityImpl();
template void Reaction<2,5>::activateReactionUnconditionalImpl();
template void Reaction<2,5>::passivateReactionImpl();
template Reaction<2,5>* Reaction<2,5>::cloneImpl(const SpeciesPtrContainerVector &spcv);

template void Reaction<4,0>::updatePropensityImpl();
template void Reaction<4,0>::activateReactionUnconditionalImpl();
template void Reaction<4,0>::passivateReactionImpl();
template Reaction<4,0>* Reaction<4,0>::cloneImpl(const SpeciesPtrContainerVector &spcv);

template void Reaction<5,2>::updatePropensityImpl();
template void Reaction<5,2>::activateReactionUnconditionalImpl();
template void Reaction<5,2>::passivateReactionImpl();
template Reaction<5,2>* Reaction<5,2>::cloneImpl(const SpeciesPtrContainerVector &spcv);

template void Reaction<4,2>::updatePropensityImpl();
template void Reaction<4,2>::activateReactionUnconditionalImpl();
template void Reaction<4,2>::passivateReactionImpl();
template Reaction<4,2>* Reaction<4,2>::cloneImpl(const SpeciesPtrContainerVector &spcv);
