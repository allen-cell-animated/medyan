//
//  Compartment.h
//  CytoSim-Experimenting
//
//  Created by Garegin Papoian on 4/21/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_Experimenting_Compartment_h
#define CytoSim_Experimenting_Compartment_h

#include <vector>
#include "Reaction.h"
#include "Species.h"
#include "SpeciesContainer.h"
#include "Composite.h"

namespace chem {
    
    class Compartment : public Composite, public SpeciesContainerVector<SpeciesDiffusing> {
    public:
        Compartment() :  Composite(), SpeciesContainerVector<SpeciesDiffusing>() {}
        virtual ~Compartment() noexcept {}
        
        virtual bool isSpeciesContainer() {return true;}
        virtual bool isReactionsContainer() {return true;}
        virtual std::string getFullName() const {return "Compartment";};
        virtual size_t countSpecies() const {return _species.size();}
        virtual size_t countReactions() const {return 0;}
        
        std::vector<SpeciesDiffusing>& species() {return _species;}
        const std::vector<SpeciesDiffusing>& species() const {return _species;}
    };



}// end of chem
//

#endif
