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
#include "Species.h"


class ParseSpecies {
private:
    std::set<Species*> _set_species;
public:
    ParseSpecies() = default;
    // SpeciesType* getSpeciesType(const std::string s){return _map_str_st.at(s);}
    Species* SpeciesProto(const std::string &name, SType stype)
    {
        SpeciesType tmpSpeciesType{name,stype};
        for (auto sptr : _set_species){
            if(sptr->getType() == tmpSpeciesType)
                return &(*sptr);
        }
        Species *s = new Species(name, stype);
        _set_species.insert(s);
        return s;
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
