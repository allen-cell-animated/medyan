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
#include <string>
#include <unordered_map>
#include <algorithm>
#include <cassert>
#include <stdexcept>

#include "common.h"
#include "utility.h"
#include "SpeciesType.h"


namespace chem {

/// We use boost::flyweights to optimize access to highly redundant data, such as SpeciesType below, 
/// which could be repeated thousands of times in many compartments. May help with cache access performance.
using namespace boost::flyweights;

class Species;
class Reaction;
class ChemSignal;

/// vr stands for vector of Reactions
typedef std::vector<Reaction*>::iterator vr_iterator; 
typedef std::vector<Reaction*>::const_iterator vr_const_iterator; 

/// vsp stands for vector of Species
typedef std::vector<Species*>::iterator vsp_iterator; 
typedef std::vector<Species*>::const_iterator vsp_const_iterator; 

    
/// Species class represents chemical molecules, tracks their copy number and can be used in [Reactions](@ref Reaction).
    
/*! This class represents chemical species, such as G-Actin. Tracks the copy number of molecules and the [Reactions](@ref Reaction)
 *   in which it is involed (@see Reaction). 
 *   @note In the future may implement the Logging and Persistence interfaces. 
 *   @note Each intantiation of Species is unique, and hence, cannot be neither copied nor moved (C++11). 
 *   This has significant implications - e.g., Species cannot be used in std::vector<Species>. Instead, 
 *   one should use either std::vector<Species*> if not owning the Species pointers, or 
 *   std::vector<std::unique_ptr<Species>> if owning the Species pointers. A special allocator can 
 *   be written such that dynamically allocated Species (through new), are arranged contigiously in memory.
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
    Species (const SpeciesType &type, species_copy_t n=0) : 
    _type(type), _n(n), _is_signaling(false) {}

    /// @param name - a string for the Species name associated with this Species. For example, "G-Actin" or "Arp2/3"
    /// @param type_enum - SType enum, such as SType::Diffusing
    /// @param n - copy number
    Species (const std::string &name, SType type_enum, species_copy_t n=0) : 
    _type(name,type_enum), _n(n), _is_signaling(false) {}

    /// deleted copy constructor - each Species is unique
    Species (const Species &r) = delete;

    /// deleted assignment operator - each Species is unique
    Species& operator=(Species&) = delete;

    /// It is required that all [Reactions](@ref Reaction) associated with this Species are destructed before this Species is destructed. 
    /// Most of the time, this will occur naturally. If not, a an assertion will ungracefully terminate the program.
    ~Species(){
        assert((_as_reactants.empty() and _as_products.empty()) && "Major bug: Species should not contain Reactions when being destroyed.");//Should not throw an exception from a destructor - that would be undefined behavior
    }
    
    /// Cloning constructs completely new Species enveloped in a unique_ptr. The copy number is set to 0.
    std::unique_ptr<Species> clone() const {return make_unique<Species>(_type,0);}

    /// Sets the copy number for this Species. 
    /// @param n should be a non-negative number, but no checking is done in run time
    /// @note The operation does not emit any signals about the copy number change.
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
    /// @param sm is the ChemSignal which will call the associated Signal (typically initiated by the 
    /// Gillespie-like simulation algorithm)
    void makeSignaling (ChemSignal &sm);
    
    /// Destroy the signal associated with this Species
    /// @param sm is the ChemSignal which manages signals
    /// @note To start signaling again, makeSignaling(...) needs to be called
    void stopSignaling (ChemSignal &sm);

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
    
    /// Print self into an iostream
    friend std::ostream& operator<<(std::ostream& os, const Species& s){
        os << s.getType() << "[" << s.getN() << "]";
        return os;
    }
};

} // end of chem namespace 
    
#endif
