//
//  SpeciesDB.h
//  CytoSim
//
//  Created by Garegin Papoian on 5/19/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_SpeciesDB_h
#define CytoSim_SpeciesDB_h

#include "Species.h"

namespace chem {

/// SpeciesDB can be used to manage the lifetime of Species. The client code will request either a Species pointer 
/// based on the existing uniq string identifier of that Species (e.g. "F:20:FA:310" - filament 20, F-actin, position 310), 
/// or a new Species can be created based on a unique identified. Cloning of existing Species is supported through this 
/// approach. 
///
/// @note It is recommended the client code uses SpeciesDB exclusively for managing Species lifetimes, or 
/// it should not be used at all, so no uncertainty arises about the Species lifetime management strategy.
/// In a large scale code, SpeciesDB should be the preferred approach.
/// @note Two SpeciesDB could be created - one for managing Species prototypes (e.g. read from the input files), and 
/// the second one the Species which will be run in a simulation.
///
/// Example:
/// @code
///    SpeciesDB sdb;
///    Species* A1 = sdb.create("A1", "A1", SType::Diffusing, 25);
///    Species* A2 = sdb.create("A2", "A2", SType::Diffusing, 25);
///    Species* A3 = sdb.create("A3", "A3", SType::Diffusing, 25);
///
///    Species *a3 = sdb1.get("A3"); // now A3 and a3 both point to the same Species
///
///    Species *c1 = sdb1.clone("C1", *A1); // A1 was cloned; @see Species::clone()
///    c1->getFullName(); // should produce "A1{Diffusing}"
///    c1->getN(); // should produce 0, because A1 was cloned.
/// @endcode
class SpeciesDB {
private:
    typedef std::unordered_map<std::string,std::unique_ptr<Species>> map_str_species;
    map_str_species _map_species; ///< a map where unique_ptr<Species> are stored, based on std::string keys
public:
    /// Defualt constructor; no Species are in the map initially
    SpeciesDB() = default;
    
    /// Copying is not allowed
    SpeciesDB(const SpeciesDB&) = delete;
    
    /// Assignment is not allowed
    SpeciesDB& operator=(SpeciesDB&) = delete;
    
    /// Create new Species, based on the last three parameters of this function. Throws std::invalid_argument if 
    /// the key was already present.
    /// @param unq_key - should be a unique string identifying this Species (supplied by the user)
    /// @param name - a string for the Species name associated with this Species. For example, "G-Actin" or "Arp2/3"
    /// @param type_enum - SType enum, such as SType::Diffusing
    /// @param n - copy number
    /// @return a pointer to the newly created Species
    Species* create(const std::string &unq_key, const std::string &name, SType type, species_copy_t N){
        auto res_pair = _map_species.insert(map_str_species::value_type( unq_key, make_unique<Species>(name, type, N)));
        if(res_pair.second==false) // insert failure -  the key was already present 
            throw std::invalid_argument("SpeciesDB::get(...) this key already exists - a bug...");
        return res_pair.first->second.get();
        //    map_species.emplace("A6",make_unique<Species>("A6", SType::Diffusing, 30)); // works with gcc 4.7
    }
    
    /// Create new Species based on the already existing Species
    /// @param proto - Species that will be cloned
    /// @param unq_key - should be a unique string identifying this Species (supplied by the user)
    /// @return a pointer to the newly created Species
    Species* clone(const std::string &unq_key, const Species &proto) {
        auto res_pair = _map_species.insert(map_str_species::value_type(unq_key, proto.clone()));
        if(res_pair.second==false) // insert failure -  the key was already present 
            throw std::invalid_argument("SpeciesDB::get(...) this key already exists - a bug...");
        return res_pair.first->second.get();
    }
    
    ///Returns a pointer to Species given a unique string key. Throws std::out_of_range on the key error.
    Species* get(const std::string &unq_key){
        auto it = _map_species.find(unq_key);
        if(it==_map_species.end())
            throw std::out_of_range("SpeciesDB::get(...) key error...");
        return it->second.get();
    }
    
    ///Return the total number of Species in the map
    size_t size() const {return _map_species.size();}
};

} // end of chem namespace 
#endif
