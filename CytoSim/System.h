//
//  System.h
//  CytoSim-Experimenting
//
//  Created by Garegin Papoian on 4/21/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_Experimenting_System_h
#define CytoSim_Experimenting_System_h

#include <vector>
#include "Reaction.h"
#include "Species.h"
#include "Compartment.h"

class System {
public:
    void setName (const std::string &name) {_name=name;}
    std::string getName () const {return _name;}
    chem::SpeciesBulk* addSpecies (const std::string name, chem::SType stype, species_copy_t n=0){
        switch(stype)
        {
            case chem::SType::Bulk:
                _species_bulk.emplace_back(new chem::SpeciesBulk(name,*this,n));
                return _species_bulk.back().get();
        }
        throw std::out_of_range("System::addSpecies: type error");
    }
private:
    std::vector<std::unique_ptr<chem::SpeciesBulk>> _species_bulk;
    std::string _name; ///< System's name
};

//namespace chem{
//    
//class SpeciesDB {
//private:
//    std::vector<Species> _vec_species;
//public:
//    size_t createSpecies(const std::string &name, SType type, species_copy_t N){
//        _vec_species.emplace_back(name,type,N);
//        return _vec_species.size();
//    }
//};

//}// end of chem namespace
#endif
