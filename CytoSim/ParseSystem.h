//
//  ParseReactions.h
//  CytoSim-Experimenting
//
//  Created by Garegin Papoian on 4/21/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_Experimenting_ParseReactions_h
#define CytoSim_Experimenting_ParseReactions_h

#include <map>
#include <stdexcept>
#include "Species.h"


class ParseSpecies {
private:
    SpeciesTypeDB _sdb;
    std::map<std::string,SpeciesType*> _map_str_st;
public:
    ParseSpecies() = default;
    SpeciesType* getSpeciesType(const std::string s){return _map_str_st.at(s);}
    bool createSpeciesTypes(const std::vector<std::string> &vs)
    {
        for(auto &s : vs){
            if(_map_str_st.count(s)>0){
                _map_str_st.clear();
                return false;
            }
            _map_str_st.insert(std::make_pair(s,_sdb.getSpeciesType(s)));
        }
        return true;
    }

    std::vector<SpeciesType*> allSpeciesType()
    {
        std::vector<SpeciesType*> ans;
        for(auto p : _map_str_st){
            ans.push_back(p.second);
        }
        std::cout << "std::vector<SpeciesType*> allSpeciesType(): " << ans.data() << std::endl;
        return ans;
    }
};

class ParseSystem {
private:
    ParseSpecies _parse_species;
public:
    ParseSystem() = default; 
    void parseSpecies(){
        _parse_species.createSpeciesTypes({"G-Actin","Arp2/3","Profilin","Motor"});
    }
    void createSpeciesTypes(const std::vector<std::string> &vs){_parse_species.createSpeciesTypes(vs);}
    SpeciesType* getSpeciesType(const std::string &s){return _parse_species.getSpeciesType(s);}
    std::vector<SpeciesType*> allSpeciesType(){return _parse_species.allSpeciesType();}
};

#endif
