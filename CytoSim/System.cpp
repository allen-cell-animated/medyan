//
//  System.cpp
//  CytoSim-Experimenting
//
//  Created by Garegin Papoian on 4/21/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//


#include "Species.h"
#include "Reaction.h"
#include "System.h"

namespace chem{
    
Species* SpeciesBulk::clone (System &other) const {
    return other.addSpecies(getName(),SType::Bulk,getN());
}

} // end of namespace chem