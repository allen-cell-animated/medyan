//
//  ReactionContainer.h
//  CytoSim
//
//  Created by Garegin Papoian on 8/31/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_ReactionContainer_h
#define CytoSim_ReactionContainer_h

#include <iostream>
#include "ReactionBase.h"
#include "Reaction.h"
    
///An abstract interface for a container of pointers to reaction objects
class ReactionPtrContainerIFace {
public:
    ///Clear the reaction container
    virtual void clear() = 0;
    
    ///Add a unique reaction pointer to this container
    virtual ReactionBase* addReactionUnique(std::unique_ptr<ReactionBase> &&Reaction) = 0;
    
    ///Remove a reaction from this container
    virtual void removeReaction(ReactionBase* Reaction) = 0;
    
    ///Find a reaction in this container
    virtual ReactionBase* findReaction (size_t index) = 0;
    
    ///Print all reactions from this container
    virtual void printReaction() {}
};

/// A concrete class implementing the ReactionPtrContainerIFace,
/// using std::vector<std::unique_ptr<ReactionBase>> as the container implementation
class ReactionPtrContainerVector : public  ReactionPtrContainerIFace {
protected:
    std::vector<std::unique_ptr<ReactionBase>> _reactions;  ///Reaction ptr container
public:
    
    ///Default constructor
    ReactionPtrContainerVector() : _reactions() {}
    
    ///Copying not allowed
    ReactionPtrContainerVector(const ReactionPtrContainerVector &) = delete;
    
    ///Assignment not allowed
    ReactionPtrContainerVector& operator=(ReactionPtrContainerVector &) = delete;

    friend void swap(ReactionPtrContainerVector& first, ReactionPtrContainerVector& second) // nothrow
    {
        // enable ADL (not necessary in our case, but good practice)
        using std::swap;
        swap(first._reactions, second._reactions);
    }

    ///Clear the container
    virtual void clear() {_reactions.clear();}

    ///Add a unique reaction ptr to this container
    virtual ReactionBase* addReactionUnique (std::unique_ptr<ReactionBase> &&Reaction) {
        _reactions.push_back(std::move(Reaction));
        return _reactions.back().get();
    }
    
    ///Add a general reaction to this container
    template<unsigned short M, unsigned short N, typename ...Args>
    ReactionBase* addReaction( Args&& ...args )
    {
        _reactions.push_back(std::unique_ptr<ReactionBase>( new Reaction<M,N>( std::forward<Args>(args)...) ));
        //        _reactions.emplace_back(make_unique(Args...));
        return _reactions.back().get();
    }
    
    ///Add a general reaction class to this container
    template<template <unsigned short M, unsigned short N> class RXN, unsigned short M, unsigned short N>
    ReactionBase* add(std::initializer_list<Species*> species, float rate)
    {
        _reactions.push_back(std::unique_ptr<ReactionBase>( new RXN<M,N>(species,rate) ));
        //        _reactions.emplace_back(make_unique(Args...));
        return _reactions.back().get();
    }
    
    ///Remove a reaction from this container
    virtual void removeReaction (ReactionBase* R) {
        auto child_iter = std::find_if(_reactions.begin(),_reactions.end(),
                                       [R](const std::unique_ptr<ReactionBase> &element)
                                       {
                                           return element.get()==R ? true : false;
                                       });
        if(child_iter!=_reactions.end())
            _reactions.erase(child_iter);
    }
    
    ///Remove all reactions that contain a certain species from this container
    virtual void removeReactions (Species* s)
    {
        for(auto &r : _reactions)
        {
            if(r->containsSpecies(s)){
                removeReaction(r.get());
                return;
            }
        }
    }
    
    ///Find a reaction by index in this container
    ///@note no check on the index
    virtual ReactionBase* findReaction (size_t index) {
        return _reactions[index].get();
    }
    
    ///Get all reactions in vector form
    std::vector<std::unique_ptr<ReactionBase>>& reactions() {return _reactions;}
    const std::vector<std::unique_ptr<ReactionBase>>& reactions() const {return _reactions;}
    
    ///Print all reactions in this container
    virtual void printReactions() {
        for(auto &r : _reactions)
            std::cout << (*r.get());
    }
    
    ///Find a similar reaction in this container (satistifes equality operator)
    ///@note returns the first similar reaction found
    virtual ReactionBase* findSimilarReaction (const ReactionBase &r) {
        auto it = std::find_if(_reactions.begin(),_reactions.end(),
                               [&r](const std::unique_ptr<ReactionBase> &element)
                               {return r==(*element);});
        if(it==_reactions.end())
            throw std::out_of_range("Reaction::findSimilarReaction(): The analogous Reaction was not found");
        return it->get();
    }
};


#endif
