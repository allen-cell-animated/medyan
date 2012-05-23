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
#include "RSpecies.h"


namespace chem {

/// We use boost::flyweights to optimize access to highly redundant data, such as SpeciesType below, 
/// which could be repeated thousands of times in many compartments. May help with cache access performance.
using namespace boost::flyweights;

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
private: //Variables
    flyweight<SpeciesType> _type; ///< This Species' type
    std::unique_ptr<RSpecies> _rspecies; 
public:
    /// Constructors 
    /// @param type - is SpeciesType associated for this Species. 
    /// @param n - copy number
    Species (const SpeciesType &type, species_copy_t n = 0) : 
    _type(type) {
        _rspecies = make_unique<RSpecies>(*this, n);
    }

    /// @param name - a string for the Species name associated with this Species. For example, "G-Actin" or "Arp2/3"
    /// @param type_enum - SType enum, such as SType::Diffusing
    /// @param n - copy number
    Species (const std::string &name, SType type_enum, species_copy_t n=0) : 
    _type(name,type_enum) {
        _rspecies = make_unique<RSpecies>(*this, n);
    }

    /// deleted copy constructor - each Species is unique
    Species (const Species &r) = delete;

    /// deleted assignment operator - each Species is unique
    Species& operator=(Species&) = delete;

    /// Virtual destructor needed for subclassing
    virtual ~Species (){}
    
    /// Return a reference to RSpecies. Notice that value copying won't be allowed 
    /// because RSpecies is not copyable.
    RSpecies& getRSpecies () {return (*_rspecies.get());}
    
    /// Return a constant reference to RSpecies. 
    const RSpecies& getRSpecies () const {return (*_rspecies.get());}
    
    /// Cloning constructs completely new Species enveloped in a unique_ptr. The copy number is copied.
    std::unique_ptr<Species> clone () const {return make_unique<Species>(_type,_rspecies->getN());}

    /// Sets the copy number for this Species. 
    /// @param n should be a non-negative number, but no checking is done in run time
    /// @note The operation does not emit any signals about the copy number change.
    void setN(species_copy_t n) {_rspecies->setN(n);}

    /// Return the current copy number of this Species
    species_copy_t getN () const {return _rspecies->getN();}
    
    /// Return SpeciesType associated with this Species
    flyweight<SpeciesType> getType () const {return _type;}

    /// Return the full name of this Species in a std::string format (e.g. "Arp2/3{Bulk}"
    std::string getFullName() const {return _type.get().getName() + "{" + _type.get().getTypeAsString() + "}";}
    
    /// Return true if this Species emits signals on copy number change
    bool isSignaling () const {return _rspecies->isSignaling();}
    
    /// Set the signaling behavior of this Species
    /// @param sm is the ChemSignal which will call the associated Signal (typically initiated by the 
    /// Gillespie-like simulation algorithm)
    void makeSignaling (ChemSignal &sm) {_rspecies->makeSignaling(sm);}
    
    /// Destroy the signal associated with this Species
    /// @param sm is the ChemSignal which manages signals
    /// @note To start signaling again, makeSignaling(...) needs to be called
    void stopSignaling (ChemSignal &sm) {_rspecies->stopSignaling(sm);}
    
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
