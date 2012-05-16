//
//  Species.h
//  CytoSim-Experimenting
//
//  Created by Garegin Papoian on 4/20/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

/**
 * @file Species.h
 * @brief this header file will contain defininition of the Species hieararchy and associated DB and helper classes.
 * @author Garegin Papoian *
 * @date 5/12/2012 
 */

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


/// @addtogroup Chemistry
/// @{

/// We use boost::flyweights to optimize access to highly redundant data, such as SpeciesType below, 
/// which could be repeated thousands of times in many compartments. May help with cache access performance.
using namespace boost::flyweights;

class SpeciesType;
class Species;
class Reaction;
class SignalingManager;

typedef std::vector<Reaction*>::iterator vr_iterator;
typedef std::vector<Species*>::iterator vsp_iterator;

enum class SType : unsigned char {
    Bulk = 0, ///< Species that have no spatial association (i.e. are "well-mixed") 
    Diffusing = 1, ///< Species that diffuse between cytosolic compartments 
    Membrane = 2, ///< Species that diffuse within a membrane 
    Filament=3, ///< Species that comprise filaments (such as F-Actin)
    Walking=4, ///< Species that can walk ("convectively") on filaments (like Myosin X)
    Motors=5 ///< Species that bind to filaments and generate forces (like Myosin II)
};

    
/// SpeciesType class represents the type of chemical Species (for example, a Species with name "Arp2/3" and SType Diffusing)

/*! This class describes the type of Species. The description is based on a name (std::string) and SType (Bulk, Diffusing, etc.). 
 *  @note In the future may implement the Logging and Persistence interfaces. 
 *  @note SpeciesType is used by Species with the flyweight pattern, since in a large system millions of unique Species may exist,
 *        but only a handful of SpeciesType objects. For example: one ("F-Actin",Filament) may correspond to millions of individual 
 *        Species.
 */
class SpeciesType {
private:
    std::string _name; ///< the descriptive name associated with this Species 
    SType _type; ///< the type of species, such as Bulk, Diffusing, etc.
    static std::vector<std::string> _vec_type_name; ///< this variable is used to help translate strings to enum values for the SType

public:
    ///Given a species name as a string and its SType, constructs the SpeciesType
    SpeciesType(const std::string &name, SType type) : _name(name), _type(type) {}

    ///Given a species name as a string and its type as a string, constructs the SpeciesType
    SpeciesType(const SpeciesType &st) : _name(st._name), _type(st._type) {}
    
    ///The move constructor. May be needed for boost::flyweight?
    SpeciesType(SpeciesType &&st) : _name(std::move(st._name)), _type(st._type) {}
    SpeciesType& operator=(const SpeciesType& st) {
        if(&st==this)
            return *this;
        _name=st._name;
        _type=st._type;
        return *this;    
    }
    
    ///Assignment operators copies all fields. May be needed for boost::flyweight?
    SpeciesType& operator=(SpeciesType&& st){
        _name=std::move(st._name);
        _type=st._type;
        return *this;
    }
    
    ///Equality operator use both species name and type to compare to SpeciesType objects
    bool operator==(const SpeciesType& species_type) const 
    {
        return species_type.getName()==_name && species_type.getType()==_type;
    }
    
    ///Returns the name associated with this SpeciesType as a string
    std::string getName() const {return _name;}

    
    ///Returns the SType associated with this SpeciesType 
    SType getType() const {return _type;}
    
    ///Returns a string representing this SpeciesType by concatanating its name and type 
    std::string getTypeAsString () const {return _vec_type_name[static_cast<int>(_type)];}
    
    ///Return if true if this SpeciesType has name and SType given as parameters to this function 
    bool is_of_type(const std::string &name, SType type) const {return name==_name && type==_type;}

    ///boost::flyweight needs a hashing function, defined here.
    friend std::size_t hash_value(SpeciesType const& species_type){
        std::size_t seed = 0;
        int type=static_cast<int>(species_type.getType());
        boost::hash_combine(seed,species_type.getName());
        boost::hash_combine(seed,type);
//        std::cout << species_type.getName()+"[" + species_type.getTypeAsString() << "] hash_value called...\n";
        return seed;
    }
};


    
/// Species class represents chemical molecules, tracks their copy number and can be used in [Reactions](@ref Reaction).
    
/*! This class represents chemical species, such as G-Actin. Tracks the copy number of molecules and the [Reactions](@ref Reaction) in which it is involed (@see Reaction). 
    @note In the future may implement the Logging and Persistence interfaces. 
 *  @note Each intantiation of Species is unique, and hence, cannot be copied. 
 */
class Species {
    /// Reactions calls addAsReactant(), removeAsReactant() - which other classes should not call
    friend class Reaction; 

private: //Variables
    std::vector<Reaction *> _as_reactants; ///< a vector of [Reactions](@ref Reaction) where this Species is a Reactant
    std::vector<Reaction *> _as_products; ///< a vector of [Reactions](@ref Reaction) where this Species is a Product
    flyweight<SpeciesType> _type; ///< This Species' type
    species_copy_t _n; ///< Current copy number of this Species
    bool _is_signaling; ///< If true, indicates a signal may be sent when a copy number of this Reaction changes
    
private:  // Methods
    /// Increases the copy number by 1. If the copy number changes from 0 to 1, calls a "callback"-like method 
    /// to activated previously passivated [Reactions](@ref Reaction), where this Species is a Reactant.
    void up() {
        _n+=1;
        if(_n==1)
            activateAssocReactions();
    }
    
    /// Decreases the copy number by 1. If the copy number changes becomes 0, calls a "callback"-like method 
    /// to passivate [Reactions](@ref Reaction), where this Species is a Reactant.
    void down() {
        _n-=1;
        if(_n == 0)
            passivateAssocReacts();
    }
    
    // \internal This methods is called by the Reaction class  during construction
    // of the Reaction where this Species is involved as a Reactant
    void addAsReactant(Reaction *r){_as_reactants.push_back(r);}
    
    // \internal This methods is called by the Reaction class during construction
    // of the Reaction where this Species is involved as a Product    
    void addAsProduct(Reaction *r){_as_products.push_back(r);}
    
    // \internal This method is called by the Reaction class during destruction
    // of the Reaction where this Species is involved as a Reactant
    void removeAsReactant(const Reaction *r) {
        auto rxit = std::find(_as_reactants.begin(),_as_reactants.end(),r);
        if(rxit!=_as_products.end()){
            _as_reactants.erase(rxit);
        }
        
    }
    // \internal This method is called by the Reaction class during destruction
    // of the Reaction where this Species is involved as a Product
    void removeAsProduct(const Reaction *r) {
        auto rxit = std::find(_as_products.begin(),_as_products.end(),r);
        if(rxit!=_as_products.end()){
            _as_products.erase(rxit);
        }
    }  
    
    /// \internal Attempts to activate previously passivated [Reactions](@ref Reaction) where this Species is involved as a 
    /// Reactant. Usually, the Reaction was first passivated, for example if the Species copy number of 
    /// one of the reactants dropeed to zero. This attempt may not succeed if there are still other
    /// reactants in the same Reaction with zero copy count.
    void activateAssocReactions();
    
    
    /// \internal Passivates all [Reactions](@ref Reaction) where this Species is among the reactants. 
    void passivateAssocReacts();

public:
    /// Constructors 
    /// @param type - is SpeciesType associated for this Species. 
    /// @param n - copy number
    Species (const SpeciesType &type, species_copy_t n=0, bool is_signaling = false) : 
    _type(type), _n(n), _is_signaling(is_signaling) {}

    /// @param name - a string for the Species name associated with this Species. For example, "G-Actin" or "Arp2/3"
    /// @param n - copy number
    Species (const std::string &name, SType type, species_copy_t n=0, bool is_signaling = false) : 
    _type(name,type), _n(n), _is_signaling(is_signaling) {}

    /// deleted copy constructor and assignment operators - each Species is unique
    Species (const Species &r) = delete;
    Species& operator=(Species&) = delete;

    /// It is required that all [Reactions](@ref Reaction) associated with this Species are destructed before this Species is destructed. 
    /// Most of the time, this will occur naturally. If not, a an assertion will ungracefully terminate the program.
    ~Species(){
        assert((_as_reactants.empty() and _as_products.empty()) && "Major bug: Species should not contain Reactions when being destroyed.");
    }
    
    /// Cloning constructs completely new Species enveloped in a unique_ptr. 
    std::unique_ptr<Species> clone() {return make_unique<Species>(_type,0);}

    /// Setters & Mutators
    /// Sets the copy number for this Species. 
    /// @param should be a non-negative number, but no checking is done in run time
    void setN(species_copy_t n) {_n=n;}

    /// Return the current copy number of this Species
    species_copy_t getN() const {return _n;}
    
    /// Return SpeciesType associated with this Species
    flyweight<SpeciesType> getType () const {return _type;}

    /// Return the full name of this Species in a std::string format (e.g. "Arp2/3{Bulk}"
    std::string getFullName() const {return _type.get().getName() + "{" + _type.get().getTypeAsString() + "}";}
    
    /// Return true if this Species emits signals on copy number change
    bool isSignaling () const {return _is_signaling;}
    
    /// Set the signaling behavior of this Species
    /// @param is the SignalingManager which will call the associated Signal (typically initiated by the 
    /// Gillespie-like simulation algorithm)
    void makeSignaling (SignalingManager &sm);
    
    /// Destroy the signal associated with this Species
    /// @param is the SignalingManager which manages signals
    /// @note To start signaling again, makeSignaling(...) needs to be called
    void stopSignaling (SignalingManager &sm);

    /// Return std::vector<Reaction *>, which contains pointers to all [Reactions](@ref Reaction) where this Species 
    /// is involved as a Reactant
    std::vector<Reaction *> ReactantReactions(){return _as_reactants;}
    
    /// Return std::vector<Reaction *>, which contains pointers to all [Reactions](@ref Reaction) where this Species 
    /// is involved as a Product
    std::vector<Reaction *> ProductReactions(){return _as_products;}
    
    /// Return std::vector<Reaction *>::iterator, which points to the beginning of all 
    /// [Reactions](@ref Reaction) where this Species is involved as a Reactant
    vr_iterator beginReactantReactions() {return _as_reactants.begin();}
    
    /// Return std::vector<Reaction *>::iterator, which points to the beginning of all 
    /// [Reactions](@ref Reaction) where this Species is involved as a Product
    vr_iterator beginProductReactions() {return _as_products.begin();}
    
    /// Return std::vector<Reaction *>::iterator, which points to the end of all 
    /// [Reactions](@ref Reaction) where this Species is involved as a Reactant    
    vr_iterator endReactantReactions() {return _as_reactants.end();}
    
    /// Return std::vector<Reaction *>::iterator, which points to the end of all 
    /// [Reactions](@ref Reaction) where this Species is involved as a Product
    vr_iterator endProductReactions() {return _as_products.end();}
    
    /// A helper function which determines whether parameter std:string name corresponds 
    /// the SpeciesType associated with this Species.
    bool is_of_species_type(const std::string &name, SType type) const {
        return _type.get().is_of_type(name,type);
    }
    
    /// Prints some helpful information about this Species. Should be used mainly for debugging 
    /// and development purposes.
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
