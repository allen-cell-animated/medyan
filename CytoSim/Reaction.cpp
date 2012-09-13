//
//  Reaction.cpp
//  CytoSim-Experimenting
//
//  Created by Garegin Papoian on 4/19/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include <boost/pool/pool_alloc.hpp>

#ifdef BOOST_MEM_POOL
    #include <boost/pool/pool.hpp>
    #include <boost/pool/pool_alloc.hpp>
#endif

#include <iostream>
#include "Reaction.h"
#include "ChemRNode.h"
#include "Composite.h"
#include "SpeciesContainer.h"

using namespace std;

namespace chem {
    
ReactionBase::ReactionBase (float rate) : 
 _rnode(nullptr), _rate(rate), _parent(nullptr)
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
    
template <unsigned short M, unsigned short N>
    void Reaction<M,N>::initializeSpecies(const std::vector<Species*> &species)
{
    assert(species.size()==(M+N) && "Reaction<M,N> Ctor: The species number does not match the template M+N");
    transform(species.begin(),species.end(),_rspecies.begin(),
              [](Species *s){return &s->getRSpecies();});

    _dependents=getAffectedReactions();
    //    cout << "Reaction::Reaction(...): " << this << endl;
    //    for (auto rr : _dependents)
    //        cout <<(*rr);
    //    cout << endl;
    //    activateReactionUnconditional();
    for(auto i=0U; i<M; ++i)
        _rspecies[i]->addAsReactant(this);
    for(auto i=M; i<(M+N); ++i)
        _rspecies[i]->addAsProduct(this);
}


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
        auto vit = std::find_if(spcv.species().cbegin(),spcv.species().cend(),
                                [molec](const std::unique_ptr<Species> &us){return us->getMolecule()==molec;});
        if(vit==spcv.species().cend())
            throw std::runtime_error("ReactionBase::Clone(): Species is not present.");
        species.push_back(vit->get());
    }
    return new Reaction<M,N>(species,_rate);
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
#endif
    
    template void Reaction<1,1>::activateReactionUnconditionalImpl();
    template void Reaction<1,1>::initializeSpecies(const std::vector<Species*> &species);
    template void Reaction<1,1>::passivateReactionImpl();
    template Reaction<1,1>* Reaction<1,1>::cloneImpl(const SpeciesPtrContainerVector &spcv);

    template void Reaction<2,1>::activateReactionUnconditionalImpl();
    template void Reaction<2,1>::initializeSpecies(const std::vector<Species*> &species);
    template void Reaction<2,1>::passivateReactionImpl();
    template Reaction<2,1>* Reaction<2,1>::cloneImpl(const SpeciesPtrContainerVector &spcv);
    
    template void Reaction<1,2>::activateReactionUnconditionalImpl();
    template void Reaction<1,2>::initializeSpecies(const std::vector<Species*> &species);
    template void Reaction<1,2>::passivateReactionImpl();
    template Reaction<1,2>* Reaction<1,2>::cloneImpl(const SpeciesPtrContainerVector &spcv);
    
    template void Reaction<2,2>::activateReactionUnconditionalImpl();
    template void Reaction<2,2>::initializeSpecies(const std::vector<Species*> &species);
    template void Reaction<2,2>::passivateReactionImpl();
    template Reaction<2,2>* Reaction<2,2>::cloneImpl(const SpeciesPtrContainerVector &spcv);
    
    
} // end of namespace
