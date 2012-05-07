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
#include <cassert>
#include <stdexcept>

#include <boost/flyweight.hpp>
#include "common.h"
#include "utility.h"


using namespace boost::flyweights;

class SpeciesType;


class Reaction;

typedef std::vector<Reaction*>::iterator VR_ITER;

enum class SType : unsigned char {Unknown=0, Bulk = 1, Diffusing = 2, Poly=3, PolyA=4, PolyM=5, ProxyA=6};

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
    std::vector<Reaction *> _as_reactants;
    std::vector<Reaction *> _as_products;
    flyweight<SpeciesType> _type;
    species_copy_t _n;
public:
    Species (const SpeciesType &type, species_copy_t n=0) : _type(type), _n(n) {}
    Species (const std::string &name, SType type, species_copy_t n=0) : _type(name,type), _n(n) {}
    Species (const Species &r) = delete;
    Species& operator=(Species&) = delete;
    ~Species(){
        assert((_as_reactants.empty() and _as_products.empty()) && "Major bug: Species should not contain Reactions when being destroyed.");
    }
    // Cloning
    Species* clone() {return new Species(_type,0);}
    // Setters & Mutators
    void setN(species_copy_t n) {_n=n;}
    void incrementN(species_copy_t delta) {_n+=delta;}
    void up() {_n+=1;}
    void down() {_n-=1;}

    void addAsReactant(Reaction *r){_as_reactants.push_back(r);}
    void addAsProduct(Reaction *r){_as_products.push_back(r);}
    
    void removeAsReactant(const Reaction *r) {
            auto rxit = std::find(_as_reactants.begin(),_as_reactants.end(),r);
            if(rxit!=_as_products.end()){
                _as_reactants.erase(rxit);
            }

    }
    void removeAsProduct(const Reaction *r) {
        auto rxit = std::find(_as_products.begin(),_as_products.end(),r);
        if(rxit!=_as_products.end()){
            _as_products.erase(rxit);
        }
        
    }    
    // Accessors 
    species_copy_t getN() const {return _n;}
    flyweight<SpeciesType> getType () const {return _type;}
    std::string getFullName() const {return _type.get().getName() + "{" + _type.get().getTypeAsString() + "}";}
    std::vector<Reaction *> ReactantReactions(){return _as_reactants;}
    std::vector<Reaction *> ProductReactions(){return _as_products;}
    VR_ITER beginReactantReactions() {return _as_reactants.begin();}
    VR_ITER beginProductReactions() {return _as_products.begin();}
    VR_ITER endReactantReactions() {return _as_reactants.end();}
    VR_ITER endProductReactions() {return _as_products.end();}
    bool is_of_species_type(const std::string &name, SType type) const {
        return _type.get().is_of_type(name,type);
    }
    //Utility
    void printSelf () const;
};

typedef std::unordered_map<std::string,std::unique_ptr<Species>> map_str_species;
    
class SpeciesContainer {
public:
    SpeciesContainer() = default;
    SpeciesContainer(const SpeciesContainer&) = delete;
    SpeciesContainer& operator=(SpeciesContainer&) = delete;
    void addSpecies(const std::string &unq_key, const std::string &name, SType type, species_copy_t N){
        auto it = _map_species.find(unq_key);
        if(it!=_map_species.end())
            throw std::invalid_argument("SpeciesContainer::getSpecies(...) this key already exists - a bug...");
        _map_species.insert( map_str_species::value_type( unq_key, make_unique<Species>(name, type, N)));
        //    map_species.emplace("A6",make_unique<Species>("A6", SType::Diffusing, 30)); // works with gcc 4.7
    }
    Species* getSpecies(const std::string &unq_key){
        auto it = _map_species.find(unq_key);
        if(it==_map_species.end())
            throw std::out_of_range("SpeciesContainer::getSpecies(...) key error...");
        return it->second.get();
    }
private:
    map_str_species _map_species;
};

    
#endif
