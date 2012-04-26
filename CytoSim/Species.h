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

enum class SType : unsigned char { Unknown=0, Bulk = 1, Diffusing = 2, Poly=3, PolyA=4, PolyM=5, ProxyA=6};

class SpeciesType{
private:
    std::string _name;
    SType _type;
    static std::vector<std::string> _vec_type_names;
public:
    //Constructors
    SpeciesType(const std::string &name, SType type) : _name(name), _type(type) {}
    SpeciesType(const SpeciesType &st) = default;
    SpeciesType(SpeciesType &&st) = default;
    SpeciesType& operator=(const SpeciesType&) = default;
    bool operator==(const SpeciesType& species_type) const {return species_type.getName()==_name && species_type.getType()==_type;}
    //Accessors
    std::string getName() const {return _name;}
    SType getType() const {return _type;}
    std::string getTypeAsString () const {return _vec_type_names[static_cast<int>(_type)];}
    //Hashing
    friend std::size_t hash_value(SpeciesType const& species_type){
        std::size_t seed = 0;
        int type=static_cast<int>(species_type.getType());
        boost::hash_combine(seed,species_type.getName());
        boost::hash_combine(seed,type);
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
    // Setters
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
    void printSelf () const;
};

class SpeciesFactory {
    
    
};

#endif

//    Species (Species &&r) :  _reactants(std::move(r._reactants)), _products(std::move(r._products)), _id(r._id), _n(r._n) {std::cout << "Species " << _id << ", Moved" << std::endl;}
