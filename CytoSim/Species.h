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
#include <unordered_map>
#include <algorithm>
#include <boost/functional/hash.hpp>
#include <boost/flyweight.hpp>
#include "common.h"

using namespace boost::flyweights;

class SpeciesType;


class ReactionBase;

enum class SType : unsigned char {Bulk = 0, Diffusing = 1, Poly=2, PolyA=3, PolyM=4, ProxyA=5};

class SpeciesType{
private:
    std::string _name;
    SType _type;
    static std::vector<std::string> _vec_type_name;
public:
    //Constructors
    SpeciesType(const std::string &name, SType type) : _name(name), _type(type) {
//        std::cout << getName()+"[" + getTypeAsString() << "] created...\n";
    }
    SpeciesType(const SpeciesType &st) : _name(st._name), _type(st._type) {
//        std::cout << getName()+"[" + getTypeAsString() << "] copied...\n";
    }
    SpeciesType(SpeciesType &&st) : _name(std::move(st._name)), _type(st._type) {        
//        std::cout << getName()+"[" + getTypeAsString() << "] moved...\n";
    }
    SpeciesType& operator=(SpeciesType&& st){
        _name=std::move(st._name);
        _type=st._type;
//        std::cout << getName()+"[" + getTypeAsString() << "] move assigned...\n";
        return *this;
    }
    SpeciesType& operator=(const SpeciesType& st) {
        _name=st._name;
        _type=st._type;
//        std::cout << getName()+"[" + getTypeAsString() << "] regular assigned...\n";
        return *this;    
    }
    bool operator==(const SpeciesType& species_type) const 
    {
//        std::cout << getName()+"[" + getTypeAsString() << "] operator==() called...\n";
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
    std::array<std::vector<ReactionBase *>,2> _reactions {};
    flyweight<SpeciesType> _type;
    species_copy_t _n=0;
public:
    Species (const SpeciesType &type) : _type(type) {}
    Species (const std::string &name, SType type) : _type(name,type) {}
    Species (const Species &r) = delete;
    Species& operator=(Species&) = delete;
    // Cloning
    Species* clone() {return new Species(_type);}
    // Setters & Mutators
    void incrementN(species_copy_incr_t delta) {_n+=delta;}
    void setN(species_copy_t n) {_n=n;}

    void addReaction(ReactionBase *r, bool left){
        if(left)
            _reactions[0].push_back(r);
        else
            _reactions[1].push_back(r);
    }
    
    void removeReaction(const ReactionBase *r) {
        for(int i=0; i<2;++i){
            auto rxit = std::find(_reactions[i].begin(),_reactions[i].end(),r);
            if(rxit!=_reactions[i].end()){
                _reactions[i].erase(rxit);
                std::cout << _type.get().getName() + "[" + _type.get().getTypeAsString() + "], Species::removeReaction, removing [" 
                << i << "] ptr->" << r << "\n";  
            }
        }
    }
    
    // Accessors 
    species_copy_t getN() const {return _n;}
    flyweight<SpeciesType> getType () const {return _type;}
    bool is_of_species_type(const std::string &name, SType type) const {
        return _type.get().is_of_type(name,type);
    }
    void printSelf () const;
};

class SpeciesFactory {
    
    
};

#endif

//    Species (Species &&r) :  _reactants(std::move(r._reactants)), _products(std::move(r._products)), _id(r._id), _n(r._n) {std::cout << "Species " << _id << ", Moved" << std::endl;}
