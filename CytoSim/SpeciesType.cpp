//
//  SpeciesType.cpp
//  CytoSim
//
//  Created by Garegin Papoian on 5/18/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include <iostream>
#include "SpeciesType.h"

namespace chem {

std::vector<std::string> vec_type_name = {"Bulk", "Diffusing", "Membrane", "Filament", "Walking", "Motors" };

std::string SpeciesType::getTypeAsString () const {
    return vec_type_name[static_cast<int>(_type)];
}
    
} // end of namespace 
