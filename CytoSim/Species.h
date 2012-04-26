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
#include "common.h"

class ReactionBase;

class SpeciesType{
public:
    SpeciesType(const std::string &type) : _type(type) {}
    SpeciesType(const SpeciesType &st) = delete;
    std::string getName() const {return _type;}
private:
    std::string _type;
};

class SpeciesTypeDB {
private: 
    std::unordered_map<std::string, SpeciesType *> _map_species_types; 
public:
    SpeciesTypeDB() = default;
    SpeciesType* getSpeciesType(const std::string &name){
        if(_map_species_types.count(name)>0) 
            return _map_species_types[name];
        SpeciesType* st = new SpeciesType(name);
        _map_species_types.insert(std::make_pair(name,st));
        return st;
    }
};

class Species {
private:
    std::array<std::vector<ReactionBase *>,2> _reactions {};
    const SpeciesType &_type;
    species_copy_t _n=0;
public:
    Species (SpeciesType &name) : _type(name) {}
    Species (const Species &r) = delete;
    Species& operator=(Species&) = delete;
    // Setters
    void incrementN(species_copy_incr_t delta) {_n+=delta;}
    void setN(species_copy_t n) {_n=n;}
    void addReaction(ReactionBase *r, bool left){
        if(left)
            _reactions[0].push_back(r);
        else
            _reactions[1].push_back(r);
    }
    void removeReaction(const ReactionBase *r){
        for(int i=0; i<2;++i){
            auto rxit = std::find(_reactions[i].begin(),_reactions[i].end(),r);
            if(rxit!=_reactions[i].end())
                _reactions[i].erase(rxit);
        }
    }
    // Accessors 
    species_copy_t getN() const {return _n;}
    void printSelf () const;
};

class SpeciesFactory {
    
    
};

#endif

//    Species (Species &&r) :  _reactants(std::move(r._reactants)), _products(std::move(r._products)), _id(r._id), _n(r._n) {std::cout << "Species " << _id << ", Moved" << std::endl;}
