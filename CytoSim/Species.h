//
//  Species.h
//  CytoSim-Experimenting
//
//  Created by Garegin Papoian on 4/20/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_Experimenting_Species_h
#define CytoSim_Experimenting_Species_h

#include <vector>
#include <array>
#include <string>
#include <algorithm>
#include <cassert>

#include <boost/flyweight.hpp>
#include "common.h"

using namespace boost::flyweights;

class SpeciesType;


class ReactionBase;

typedef std::vector<ReactionBase*>::iterator VRB_ITERATOR;

enum class SType : unsigned char {Bulk = 0, Diffusing = 1, Poly=2, PolyA=3, PolyM=4, ProxyA=5};

class SpeciesType{
private:
    std::string _name;
    SType _type;
    static std::vector<std::string> _vec_type_name;
public:
    //Constructors
    SpeciesType(const std::string &name, SType type) : _name(name), _type(type) {}
    SpeciesType(const SpeciesType &st) : _name(st._name), _type(st._type) {}
    SpeciesType(SpeciesType &&st) : _name(std::move(st._name)), _type(st._type) {}
    SpeciesType& operator=(const SpeciesType& st) {
        if(&st==this)
            return *this;
        _name=st._name;
        _type=st._type;
        return *this;    
    }
    SpeciesType& operator=(SpeciesType&& st){
        _name=std::move(st._name);
        _type=st._type;
        return *this;
    }
    bool operator==(const SpeciesType& species_type) const 
    {
        return species_type.getName()==_name && species_type.getType()==_type;
    }
    //Accessors
    std::string getName() const {return _name;}
    SType getType() const {return _type;}
    std::string getTypeAsString () const {return _vec_type_name[static_cast<int>(_type)];}
    bool is_of_type(const std::string &name, SType type) const {return name==_name && type==_type;}
    //Hashing
    friend std::size_t hash_value(SpeciesType const& species_type){
        std::size_t seed = 0;
        int type=static_cast<int>(species_type.getType());
        boost::hash_combine(seed,species_type.getName());
        boost::hash_combine(seed,type);
//        std::cout << species_type.getName()+"[" + species_type.getTypeAsString() << "] hash_value called...\n";
        return seed;
    }
};


class Species {
private:
    std::vector<ReactionBase *> _freactions;
    std::vector<ReactionBase *> _breactions;
    flyweight<SpeciesType> _type;
    species_copy_t _n;
public:
    Species (const SpeciesType &type, species_copy_t n=0) : _type(type), _n(n) {}
    Species (const std::string &name, SType type, species_copy_t n=0) : _type(name,type), _n(n) {}
    Species (const Species &r) = delete;
    Species& operator=(Species&) = delete;
    ~Species(){
        assert((_freactions.empty() and _breactions.empty()) && "Major bug: Species should not contain Reactions when being destroyed.");
    }
    // Cloning
    Species* clone() {return new Species(_type,0);}
    // Setters & Mutators
    void setN(species_copy_t n) {_n=n;}
    void incrementN(species_copy_t delta) {_n+=delta;}
    void up() {_n+=1;}
    void down() {_n-=1;}

    void addFReaction(ReactionBase *r){_freactions.push_back(r);}
    void addBReaction(ReactionBase *r){_breactions.push_back(r);}
    
    void removeFReaction(const ReactionBase *r) {
            auto rxit = std::find(_freactions.begin(),_freactions.end(),r);
            if(rxit!=_breactions.end()){
                _freactions.erase(rxit);
            }

    }
    void removeBReaction(const ReactionBase *r) {
        auto rxit = std::find(_breactions.begin(),_breactions.end(),r);
        if(rxit!=_breactions.end()){
            _breactions.erase(rxit);
        }
        
    }

    
    // Accessors 
    species_copy_t getN() const {return _n;}
    std::vector<ReactionBase *> getFReactions(){return _freactions;}
    VRB_ITERATOR beginFReactions() {return _freactions.begin();}
    VRB_ITERATOR beginBReactions() {return _breactions.begin();}
    VRB_ITERATOR endFReactions() {return _freactions.end();}
    VRB_ITERATOR endBReactions() {return _breactions.end();}
    flyweight<SpeciesType> getType () const {return _type;}
    bool is_of_species_type(const std::string &name, SType type) const {
        return _type.get().is_of_type(name,type);
    }
    std::string getFullName() const {return _type.get().getName() + "[" + _type.get().getTypeAsString() + "]";}
    void printSelf () const;
};

class SpeciesFactory {
    
    
};

    
#endif
