//
//  ParseReactions.h
//  CytoSim-Experimenting
//
//  Created by Garegin Papoian on 4/21/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_Experimenting_ParseReactions_h
#define CytoSim_Experimenting_ParseReactions_h

#include <set>
#include <map>
#include <stdexcept>
#include <boost/bind.hpp>
#include "Species.h"

//struct has_species_type {
//     bool operator()(Species *x,  const SpeciesType &y) {return x->getType() == y;}
//}; 

class ParseSpecies {
private:
    std::set<Species*> _set_species;
public:
    ParseSpecies() = default;
    // Given a Species name and type, look for already existing similar Species in the internal database. If found
    // return that Species, otherwise create a new one and return
    Species* SpeciesProto(const std::string &name, SType stype)
    {
        for (auto &sptr : _set_species){
            if(sptr->is_of_species_type(name,stype))
                return &(*sptr);
        }
        
        // If such Species does not exist, create a new one
        Species *s = new Species(name,stype);
        _set_species.insert(s);
        return s;

//       ##### The commented out codes below show C++11 and pre-C++11 ways of accomplishing the same as above.
        
//       ##### First approach is based on using a lambda function and boost::bind:
//        
//        auto sptr1 = std::find_if(_set_species.begin(),_set_species.end(), 
//                                  [&tmpSpeciesType](Species *x)
//        {
//            return x->getType() == tmpSpeciesType;
//        });
//
//        if(sptr!=_set_species.end())
//            return *sptr;
        

//        ##### Second alternative approach is based on an explicit predicate (i.e. without C++11 lambda)
//
//        auto cmp_pred = boost::bind<bool>(has_species_type(),_1,boost::cref(tmpSpeciesType));
//        auto sptr = std::find_if(_set_species.begin(),_set_species.end(), cmp_pred);
        
//        Second alternative approach to do the same        
        
    }
};

class ParseSystem {
private:
    ParseSpecies _parse_species;
public:
    ParseSystem() = default; 
    void parseSpecies(){
        _parse_species.SpeciesProto("G-Actin",SType::Diffusing);
    }
    Species* SpeciesProto(const std::string &name, SType stype){return _parse_species.SpeciesProto(name, stype);}
};

#endif
