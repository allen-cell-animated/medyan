//
//  Protofilament.h
//  CytoSim
//
//  Created by Garegin Papoian on 6/2/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_Protofilament_h
#define CytoSim_Protofilament_h

#include "Composite.h"
#include "Species.h"

namespace chem {
    
class ProtoFilament : public Composite {
private:
    SpeciesDiffusing _proto_species;
public:
    ProtoFilament(const std::string &species_name, size_t len = 0) : Composite(), _proto_species(species_name) 
    {
        for(size_t i=0; i<len; ++i){
            auto species = std::unique_ptr<SpeciesDiffusing>(new SpeciesDiffusing(_proto_species));
            addSpeciesUnique(std::move(species));
            this->species().back()->setParent(this);
        }
    }

    virtual void grow(size_t len) {
        for(size_t i=0; i<len; ++i){
            auto species = std::unique_ptr<SpeciesDiffusing>(new SpeciesDiffusing(_proto_species));
            addSpeciesUnique(std::move(species));
            this->species().back()->setParent(this);
        }
    }
    
    virtual std::string getFullName() const {
        return "ProtoFilament:[" + _proto_species.getName() + "]";
    }

};

} // end of chem


#endif
